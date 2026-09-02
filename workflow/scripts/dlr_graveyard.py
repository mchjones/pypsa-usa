# old code from in-flow dlr calculations
## if you're doing raw DLR calculations with clustered bus locations in the PyPSA-USA workflow
constant_inputs = {
    "T_s": 100,
    "T_avg": 100,
    "em": 0.8,
    "alpha": 0.8
}
static_rating = { # following the ieee guidance for selection of weather variables for base ratings: E.1 Base Ratings from the IEEE standard and 4.6.1.1 [1]
    "v_w": 0.61, 
    "phi": 90, 
    "Q_s": 1000,
    "T_a": 40,
    "hour": 11, 
    "day": 161
}
solar_coeff_data = {
    0: -42.2391,
    1: 63.8044,
    2: -1.9220,
    3: 3.46921e-2,
    4: -3.61118e-4,
    5: 1.94318e-6,
    6: -4.07608e-9
}

def load_raw_WUS_data(planning_horizon: int) -> xr.Dataset:
    #Loads WUS data for given planning horizon
    base_path = snakemake.config['cf_path']
    logger.info(f"Loading raw WUS data for planning horizon {planning_horizon} from {base_path}...")
    file_path = Path(base_path) / f"{planning_horizon}/regrid_{planning_horizon}_ssp370_d02.nc"
    
    if not file_path.exists():
        raise FileNotFoundError(f"WUS data file not found at: {file_path}")
        
    return xr.open_dataset(file_path)

def calculate_dlr(n):
    if snakemake.wildcards.dlr != "none":
        logger.info("Starting dlr calculation...")

        # get line and bus data
        elevation = pd.read_csv(snakemake.input.elev_ref,header=0,index_col=0)
        line_data = get_line_info(n,elevation)
        bus_data = get_bus_data(n)
        lines = n.lines.index.values
        num_lines = len(lines)

        line_inputs, bus_inputs, bus_list = get_line_dlr_inputs(line_data)

        baseline = calc_system_baseline(line_inputs,static_rating,constant_inputs,solar_coeff_data)

        # set assumptions - eventually put this into config maybe
        wind_fix = False # if True, assumes v = 0.61 m/s, phi = 90 deg (based on IEEE standard)
        phi_fix = True # if True, uses WUS to calculate magnitude, but fixes phi = 0 deg

        # for each planning horizon, load wus data and calc dlrs
        year_relative = {}
        for i, year in enumerate(snakemake.config['scenario']['planning_horizons']):
            wus_data_horizon = load_raw_WUS_data(year)
            wus_data_prev = load_raw_WUS_data(int(year)-1)
            wus = get_fixed_year(wus_data_horizon,wus_data_prev,year)
            snap = gen_year_no_leap(year)

            wus_extracted = wus.sel(
                lat=xr.DataArray(bus_data['lat'], dims='points'),
                lon=xr.DataArray(bus_data['lon'], dims='points'),
                method='nearest'
            )[['T2', 'U10', 'V10', 'SWDNB']]
            wus_extracted = wus_extracted.assign_coords(Times=np.arange(0,8760,1))

            year_dlr = get_line_dlrs_eff(wus_extracted,constant_inputs,line_inputs,solar_coeff_data,bus_inputs,bus_list,num_lines,wind_fix,phi_fix)
            rel = get_relative(year_dlr,baseline)
            rel.index = pd.DatetimeIndex(snap)
            year_relative[year] = rel

        logger.info("Yearly DLR calculations complete, stacking")
        all_relative = pd.concat([year_relative[year] for year in sorted(year_relative.keys())])

        if snakemake.wildcards.dlr == "dlr": # truncate at 1.3
            dlr = trunk(all_relative,1.3)
            dlr = dlr.round(3)
            dlr.to_csv(snakemake.output.dlr_path)
            logger.info(f"{snakemake.wildcards.dlr} exported to {snakemake.output.dlr_path}.")
        elif snakemake.wildcards.dlr == "derate": # truncate at 1
            derate = trunk(all_relative,1)
            derate = derate.round(3)
            derate.to_csv(snakemake.output.dlr_path)
            logger.info(f"{snakemake.wildcards.dlr} exported to {snakemake.output.dlr_path}.")
        elif snakemake.wildcards.dlr == "slr": # all 1's
            slr = pd.DataFrame(np.ones_like(all_relative), 
                      index=all_relative.index, 
                      columns=all_relative.columns)
            slr.to_csv(snakemake.output.dlr_path)
            logger.info(f"{snakemake.wildcards.dlr} exported to {snakemake.output.dlr_path}.")
           
    else:
        logger.info(f"DLR is set to {snakemake.wildcards.dlr}. No DLR added")
## raw calculations direct in PyPSA-USA workflow
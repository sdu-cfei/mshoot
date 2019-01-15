Learning period:    2017-03-21 00:00:00 - 2017-04-05 00:00:00  (15 days)
Validation period:  2017-04-05 00:00:00 - 2017-04-21 00:00:00  (16 days)


File:              | Description:
===================|============================================================
areas.txt          | Room areas [m2]
volumes.txt        | Room volumes [m2]
radiators.txt      | Room radiator capacities [W]
df_t.csv           | Room temperature [degC]
df_blinds.csv      | Blind position [%] (100 = up, 0 = down)
df_co2.csv         | Room CO2 [ppm]
df_FTI1.csv        | AHU sensor, temp. after heating coil (see AHU.png) [degC]
df_FTI2.csv        | AHU sensor, temp. before passing though HX (see AHU.png) [degC] (can be used to validate outdoor temperature)
df_hpos.csv        | Radiator valve position [byte 0-255]
df_occ.csv         | Camera-based occupancy [persons]
df_pir.csv         | Motion sensor [binary 0-1]
df_solrad.csv      | Solar radiation [W/m2]
df_tamb.csv        | Outdoor temperature [degC] (note that the sensor is likely not shaded, compare with FTI2)
df_windspeed.csv   | Wind speed [m/s]
df_vav.csv         | Damper position (per room) [%]
df_el_light.csv    | Lighting el. consumpt. (sum of 3 light rows) [kWh] - only e22-511-2 - likely incomplete!
df_plug_loads.csv  | Electricity consumption (plug loads only) [kWh]    - only e22-511-2
df_heat_supply.csv | Radiator heat supply [kWh]                         - only e22-511-2


Black-box model
===============
Suggested inputs:               hpos, occ, solrad, tamb, vav
Possible additional inputs:     FTI1, blinds
Obligatory output:              t

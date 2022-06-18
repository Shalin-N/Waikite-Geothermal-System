# Project Background
You have been retained by Tūhourangi Ngāti Wāhiao to undertake a computer modelling study of the Rotorua Geothermal Field and recovery of the Waikite Geyser. To support your study, the following data are available:
- Total extraction rate data from the Rotorua geothermal system, and a partition corresponding to a Rhyolite formation that has seen continuing extraction post borehole closure.
- Pressure and temperature monitoring data from near Whakarewarewa.
- Reservoir engineering reports which indicate that the supply of hot water (rate and temperature) into the geysers has been disrupted but is recovering. The supply of cold water is unchanged.

## Project expectations:
You should undertake a computer modelling study that will assist decision-making during the resource consent hearing, in particular addressing the noted concerns of other stakeholders where they are relevant to the study. The model you develop should be defensible, reflective of reality, and take appropriate account of uncertainty. You will be required to communicate the model findings in both oral and written formats.

Following list of files can be found in the ***code*** folder. All diagrams will be produced when the ***main.py*** script is ran. Every other folders are irrelevant for the assessment of our project. ***pressure.py*** and ***temperature.py*** python scripts are called within ***main.py*** python script and are used for calculations. ***tests.py*** is used for unit testing. 

## The list of diagrams and their description:

### Given

- **given_data_1.png** : plot of given data (water level and total production rate, q1 with _reinjection_ considered).

- **given_data_2.png** : plot of total production rate, q1 with _reinjection_ considered, and water level converted to pressure (in units of bar).

- **given_data_3.png** : plot of given data (temperature and total production rate, q1 with _reinjection_ considered).

### Assume 

- **approximate_area.png** : the reservoir area estimation. 

- **Pressure_sketch.png** : conceptual diagram of reservoir (pressure).

- **Temperature_sketch.png** : conceptual diagram of reservoir (temperature). 

### Working

- **pressure_validation.png** : plot of pressure benchmarking and pressure convergence analysis.

- **temperature_validation.png** : plot of temperature benchmarking and temperature convergence analysis.

### Suitable

- **first_fit.png** : plot of initial pressure best fit model (cp = 0, which disregards slow drainage).

- **first_pressure_misfit.png** : plot of initial pressure best fit model and misfit to data quantified.

- **first_temperature_misfit.png** : plot of initial temperature best fit model and misfit to data quantified.

### Improve

- **best_fit.png** : plot of best fit pressure and temperature LPM ODE models.

- **best_pressure_misfit.png** : plot of best fit pressure LPM ODE model and misfit to data quantified. 

- **best_temperature_misfit.png** : plot of best fit temperature LPM ODE model and misfit to data quantified. 

### Use 

- **pressure_forecast.png** : plot of pressure forecast for different production rate scenarios.

- **pressure_forecast_supplement.png** : plot of rate of pressure change forecast for different production rate scenarios. The dashed line shows the rate of pressure change in which the Rotorua geothermal system is not recovering (rate of pressure change must be greater than 0 at all times to show recovery). 

- **temperature_forecast.png** : plot of temperature forecast for different production rate scenarios (temperature is not directly dependent on production rates, but it is dependent on pressure values that is dependent on production rates).

- **temperature_forecast_supplement.png** : plot of rate of temperature change forecast for different production rate scenarios. The dashed line shows the rate of temperature change in which the Rotorua geothermal system is not recovering (rate of temperature change must be greater than 0 at all times to show recovery).

### Unknown

- **pressure_forecast_uncertainty.png** : plot of pressure forecast with uncertainty.

- **pressure_forecast_uncertainty_supplement.png** : plot of rate of change of pressure forcast with uncertainty.

- **temperature_forecast_uncertainty.png** : plot of temperature forecast with uncertainty.

- **temperature_forecast_uncertainty_supplement.png** : plot of rate of change of temperature forecast with uncertainty.


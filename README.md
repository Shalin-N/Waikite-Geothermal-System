# Project Background
You have been retained by Tūhourangi Ngāti Wāhiao to undertake a computer modelling study of the Rotorua Geothermal Field and recovery of the Waikite Geyser. To support your study, the following data are available:
- Total extraction rate data from the Rotorua geothermal system, and a partition corresponding to a Rhyolite formation that has seen continuing extraction post borehole closure.
- Pressure and temperature monitoring data from near Whakarewarewa.
- Reservoir engineering reports which indicate that the supply of hot water (rate and temperature) into the geysers has been disrupted but is recovering. The supply of cold water is unchanged.

## Project expectations:
You should undertake a computer modelling study that will assist decision-making during the resource consent hearing, in particular addressing the noted concerns of other stakeholders where they are relevant to the study. The model you develop should be defensible, reflective of reality, and take appropriate account of uncertainty. You will be required to communicate the model findings in both oral and written formats.

## Code
Following list of files can be found in the ***code*** folder. All diagrams will be produced when the ***main.py*** script is ran. Every other folders are irrelevant for the assessment of our project. ***pressure.py*** and ***temperature.py*** python scripts are called within ***main.py*** python script and are used for calculations. ***tests.py*** is used for unit testing. 

# 1. Why?
Rotorua city was built atop the Rotorua geothermal system, and for decades its residents drilled shallow bores to use free hot water for heating and bathing. Over time though, the exploitation of this field caused:
- Pressure and groundwater level to decline
-	Geysers supported by this system  to stop erupting, like the Waikite Geyser at Whakarewarewa

The Bay of Plenty Regional Council has received an application from the Rotorua City Council, which wishes to extend the borehole moratorium. There is currently limited production from existing bores and a borehole exclusion zone around Whakarewarewa. Whether the council approves or declines this application depends upon the submissions of the stakeholders:
-	Rotorua City Council: (RCC) 
-	The local iwi (Tūhourangi Ngāti Wāhiao)
-	The local Chamber of Commerce: representing local hotels


The RCC has proposed to continue the moratorium in its current form. The Tūhourangi Ngāti Wāhiao oppose this as they wish to see all boreholes closed and the local chamber of commerce representing local hotels who want to install geothermally heated baths. Based on the interest of stakeholders, possible outcomes are:
-	The application is accepted, and the moratorium is extended in its current form.
-	The application is rejected, and the borehole operation must cease.
-	The application is rejected, and the moratorium is lifted, allowing for additional boreholes.

# 2.	How?
A computer modelling study can be used to understand how changes in operation (borehole extraction and cold water reinjection) affects the pressure and temperature, and surface features of a geothermal system. These insights could be used to model the economics (tourism) of the region, but that will be considered out of the scope of this study. 

Tūhourangi Ngāti Wāhiao has retained us to perform a computer modelling study of the Rotorua Geothermal field and, subsequently, the Waikite Geyser. Our model aims to simulate how pressure and temperature in the system change in response to different levels of operation (borehole extraction and cold water reinjection). This model will predict reservoir pressure, temperature, and surface features for different operation rates.

The results from this modelling study will provide insight into:
-	what level of operation can be sustained without adversely impacting geothermal surface features
-	How the system will react to the addition/removal of borewells
-	How the system will respond to varying rates of operation
-	How the system will recover if left in its current state

Previous studies suggest that surface features are negatively impacted when the reservoir pressure drops below 0.02 bar.

# 3. Given?
The Rotorua Geothermal System can be modelled as a shallow or deep recharge model. Beneath Rotorua lake, there exist columns that allow cold water (10C) at constant hydrostatic pressure to flow into the geothermal system.

A lumped parameter model can be used to conceptualise a geothermal reservoir as a single block. The pressure changes within the system in response to fluid extraction and inflow of adjacent groundwater (Fradkin, LJ, Storey, & McNabb, 1981). The temperature changes within the system according to cold water inflow and conduction. These can be represented as an ODE and are coupled together using pressure. 

![image](https://user-images.githubusercontent.com/85419997/174466096-8538ba9d-cc14-4712-bfc7-052f97e1d4c9.png)
<br><i>Figure 1: Averge water level data from monitoring wells from near Whakarewarewa (blue) plotted against extraction rate (black) from 1950 to 2014</i>

![image](https://user-images.githubusercontent.com/85419997/174466111-4c234e5c-a32c-4224-8be4-f04944fa7737.png)
<br><i>Figure 2: Averge pressure data from monitoring wells from near Whakarewarewa (blue) plotted against extraction rate (black) from 1950 to 2014</i>

![image](https://user-images.githubusercontent.com/85419997/174466122-ca0d0331-0bb2-4962-9b54-fc1ad3d8c217.png)
<br><i>Figure 3: temperature readings reservoir engineering reports (red) plotted against extraction rate (black) from 1950 to 2014</i>

 We have been given data from a monitoring well from near Whakarewarewa. The data shows that from 1984 when the moratorium was imposed, the water level slowly began to recover. This water level can be converted to pressure readings representing the geothermal reservoir. We have also been given reservoir engineering reports which indicate that the supply of hot water (rate and temperature) into the geysers had been disrupted but is recovering. The supply of cold water is unchanged. We can see the decline in water level and temperature coincide with the increase in extraction rate and the recovery with the decrease in extraction rate.

Other useful data such as water level readings from 1950 to 1984, Annual rainfall data from 1950 to 2014 and the approximate size of the geothermal reservoir is 15 to 28 km2 have been discovered in the literature.

#### References
<br>Allis, R. G., & Lumb, J. T. (1992). The Rotorua geothermal field, New Zealand: its physical setting, hydrology, and responce to exploitation. Geothermics 21, 7-24.
Bradford, E. (1992). Pressure changes in Rotorua geothermal aquifers, 1982-90. Geothermics 21, 231-248.

Cody, A. D., & Lumb, J. T. (1992). Changes in thermal activity in the Rotorua geothermal field. Geothermics 21, 215-230.

Fradkin, LJ, Storey, M., & McNabb, A. (1981). On Identification and Validation of Some Geothermal Models. Water Resources Research 17, 929-936.

Ratouis, T. M., O’Sullivan, M. J., Alcaraz, S. A., & O’Sullivan, J. P. (2017). The effects of seasonal variations in rainfall and production on the aquifer and surface features of Rotorua geothermal field. Geothermics 69, 165-188.

Rinehart, J. S. (1980). Geysers and Geothermal Energy. New York: Springer-Verlag.

Saptadji, N., O’Sullivan, J., Krzyzosiak, W., & O’Sullivana, M. (2016). Numerical modelling of Pohutu geyser, Rotorua, New Zealand. Geothermics 64, 401-409.

Scott, B. J., Gordon, D. A., & Cody, A. D. (2005). Recovery of Rotorua geothermal field, New Zealand: Progress, issues and consequences. Geothermics 34, 161-185.

Scott, B. J., Mroczek, E. K., Burnell, J. G., Zarroukc, S. J., Seward, A., Robson, B., & Graham, D. J. (2016). The Rotorua Geothermal Field: An experiment in environmental management. Geothermics 59, 294-310.

# 4. Assume?
The physics relevant to a geothermal system are:
- Conservations of mass, which governs changes in fluid pressure. Fluid flow is by diffusion and is described by Darcy’s Law 
-	Conservation of energy, which governs changes in fluid and rock temperature. Heat transfer is both advective (heat carried by moving fluid) and conductive (heat diffusing through rock). Conductive heat is described by Fourier’s Law
-	Conservation of momentum, which governs changes in stress and deformation. Most rocks are well-described by linear elasticity

Our concern here is the understanding of fluid pressure and temperature in the Rotorua geothermal System. Pressure and temperature affect whether a given extraction rate is sustainable long term and whether there are likely to adverse impacts on surface features such as the Waikite Geyser.

For this model we will consider the conservation of mass and energy. We shall use a LPM like that described by (Fradkin, LJ, Storey, & McNabb, 1981) to model the pressure change and a LPM obtained by rewriting the conservation of thermal energy equation to model the temperature change. Our 2 models will be solved in separate solution spaces and be coupled using pressure predictions from the pressure LPM.

Choosing domain:
The LPM approach requires that we represent the entire geothermal reservoir as a large single control volume. We assumed the Rotorua Geothermal System is 15.8 km^2. The results of the model are to be applied only within this domain. The time-domain shall be from the start of the decline of the health of the Rotorua Geothermal System (1950) until the end of our data (2014).

To develop an accurate LPM, we will need to constrain the additional conditions:
-	The initial pressure of Rotorua Geothermal System in 1950
-	The initial temperature of Rotorua Geothermal System in 1950
-	Pressure values at the recharge source
-	Temperature value at the recharge source
-	Temperature value of cold-water inflow

The effect of adjacent groundwater to the reservoir is included as a term in the LPM ODE. This is a boundary condition. The strength of this boundary condition is unknown and will be inferred by calibration to the data. 

![image](https://user-images.githubusercontent.com/85419997/174466188-1c35d6a4-2305-4199-986c-76b059edceb5.png)
<br><i>Figure 3: sketch representing the pressure of Rotorua geothermal reservoir</i>

![image](https://user-images.githubusercontent.com/85419997/174466201-321e2ad6-5678-4fb6-9183-cddd6dbcc5dc.png)
<br><i>Figure 4: sketch representing the temperature of Rotorua geothermal reservoir</i>

Parameters:
Specified in formulate

# 5. Formulate?
![image](https://user-images.githubusercontent.com/85419997/174466290-2f7e002d-c966-49ee-bd32-ffb68bab0385.png)

![image](https://user-images.githubusercontent.com/85419997/174466301-9a907b8e-12e5-4b58-9e4f-111492820bd3.png)

# 6. Working?
We have checked that the Improved Euler implementation is working correctly by running a unit test that compares the output to a by-hand example. The code is given in test.py.

We have benchmarked the numerical solution of the LPM against the analytical solution obtained for constant production rate. The results are shown below and indicate the solver is working as expected. For the temperature numerical solution, the rate of change is quicker, but they both stabilise at the same point. We also see evidence of instability for large time steps in both temperature and pressure. We can avoid this error by using a time step of 1.

![image](https://user-images.githubusercontent.com/85419997/174466323-c87c5cc4-b9b3-4a47-af0f-e9f2f0ffa8fb.png)
<br><i>Figure 5: comparison of numerical solution to analytical (left) and plot of instability for time steps (right) fore pressure</i>

![image](https://user-images.githubusercontent.com/85419997/174466331-92bd388d-cb34-46d4-879b-044832b8c7d3.png)
<br><i>Figure 6: comparison of numerical solution to analytical (left) and plot of instability for time steps (right) fore pressure</i>

# 7. Suitable?
I have calibrated the LPM model without slow drainage to the given data using curve_fit, which uses non-linear least squares to fit a function. The best fit model has pressure parameters a=110.366,   b=1,c=0 and temperature parameters a=4.189e-05,   b=0.086

![image](https://user-images.githubusercontent.com/85419997/174466352-77c930de-8c64-4e82-947e-f289b165d643.png)
<br><i>Figure 7: Best-fit model for pressure without slow drainage and unknown data from 1950-1984 (right) compared against pressure data (blue). Plot of model misfit against data (left)</i>

![image](https://user-images.githubusercontent.com/85419997/174466368-b650457a-60f0-4318-8e79-9b8762b6e0fc.png)
<br><i>Figure 8: Best-fit model for temperature (right) compared against temperature data (red). Plot of model misfit against data (left)</i>

In regard to pressure:
This model does a good job of fitting the final pressure data but does not accurately capture the rise and rebound point. The misfit is correlated in time and indicates the model’s departure from reality. The period from 1984 to 1995 has an extreme misfit of 0.2 but reduces to < 0.05 as time progresses. This may be due to the absence of slow drainage physics and the uncertainty in data from 1950 to 1984.

In regard to temperature:
This model does a good job fitting the trend of temperature change and captures the increasing pressure, albeit a bit lower than the given data. The misfit is correlated in time and indicates the model’s departure from reality. This model has an average misfit of 1.5 °C. This misfit may be because we only have a small data set from 1984 to 2014 with measurements every five years. 

# 8. Improve?
Mass extraction and groundwater recharge are sufficient physics to describe the major trends in the data, but for the more subtle trends, some changes that may improve the model are:
-	Inclusion of a slow drainage term in the LPM. 
-	Inclusion of data from before 1984.
-	Inclusion of reinjection of cold water

I have calibrated the model a second time, now allowing the parameter 𝑐 to be variable approximating the data from before 1984 according to an arc of a circle (This has been omitted from the misfit graph because it is not accurate data). The best fit model now has a=8.688,   b=0.105,c=66.682  and a=5.104e-07,   b=0.084.

![image](https://user-images.githubusercontent.com/85419997/174466413-441a2081-f5bc-4ff7-95c1-c32b357ad7f4.png)
<br><i>Figure 9: Best-fit model for pressure with slow drainage and unknown data from 1950-1984 (right) compared against pressure data (blue). Plot of model misfit against data (left)</i>

![image](https://user-images.githubusercontent.com/85419997/174466424-f7d09930-0d49-402a-99fb-ad68c5e0c162.png)
<br><i>Figure 10: Best-fit model for temperature (right) compared against temperature data (red). Plot of model misfit against data (left)</i>

In regard to pressure:
This model does a better job fitting all pressure data points and accurately models the stabilisation of pressure. The misfit is correlated in time and indicates the model’s departure from reality. The misfit has now been reduced to 0.05 bar at its extreme but < 0.03 bar on average. This may be because our model approximates the Rotorua Geothermal System as a single CV and cannot capture finer detail such as separate cavities. 

In regard to temperature:
This model does a similar job of fitting the trend of temperature change and captures the increasing pressure. The misfit is correlated in time and indicates the model’s departure from reality. This model has an average misfit of 1.5 °C, similar to the previous model. This misfit may be due to how seasons affect the temperature of cold-water inflow and conduction, which our model does not account for.

# 9.	Use?
We will use the model calibrated to parameters from section 8 to consider four “what-if” scenarios that will end in 2050.

These are:
<br>i.	Maintain extraction at the current rate.
<br>ii.	Stop all extraction.
<br>iii.	Doubled extraction of the current rate.
<br>iv.	Half extraction of the current rate.

These scenarios represent the possible outcomes of the application and identify the range of change in the Rotorua Geothermal System temperature and pressure and hence a possible impact on the Waikite Geyser. We have modelled the rate of change to be non-constant for a more realistic prediction as it is unlikely that all operations could be halted or doubled instantly.

![image](https://user-images.githubusercontent.com/85419997/174466456-01f33eb5-aab6-4f3d-99a9-d58de2282911.png)
<br><i>Figure 11: pressure model predictions for four scenatio plotted with given pressure data (blue)</i>

![image](https://user-images.githubusercontent.com/85419997/174466467-8c5e073e-e7d4-4418-b4fe-6f46318afc13.png)
<br><i>Figure 12: pressure model predictions for four scenatio plotted with given pressure data (blue)</i>

![image](https://user-images.githubusercontent.com/85419997/174466507-460c53b5-c8ad-4b0d-a6d2-a15cf782f284.png)
<br><i>Figure 13: pressure model rate of change of predictions for the four scenarios</i>

![image](https://user-images.githubusercontent.com/85419997/174466515-0fef4ffe-acbe-4193-a884-1c4e1652e87a.png)
<br><i>Figure 14: pressure model rate of change of predictions for the four scenarios</i>

The final pressures for each scenario are:	            
i.	0.0288 bar				   	
<br>ii.	0.0478 bar				   	
<br>iii.	0.0095 bar				  	
<br>iv.	0.0384 bar				  	

The final temperatures for each scenario are:
<br>i. 	145.11 °C
<br>ii.	146.13 °C
<br>iii.	144.095 °C
<br>iv.	145.62 °C

In regards to pressure:
For all scenarios, the rate of change of pressure reaches a stable state. The stop & half current extraction scenarios see an increase in pressure. The maintained extraction scenario reaches a stable state. The double extraction scenario steadily declines. 

In regards to temperature:
For all scenarios, the rate of change of temperature declines to a stable state scaled by the extraction rate. The stop & half current extraction scenarios see a steady increase in temperature. The maintained extraction scenario reaches a stable state. The double operation scenario initially increases but then begins to decline.

# 10. Unknown?
There is uncertainty in the temperature and pressure data due to the missing data from 1950-1984, the variation in annual rainfall and temperatures associated with seasonal change. We approximated that the water level measurements had a standard deviation equal to yearly rainfall, and the temperature measurements had a standard deviation of 2 °C due to seasons temperature fluctuation. By minimising the sum of least squares of an error where sigma is equal to the standard deviations mentioned above, we can reasonably estimate the model parameters.

We used the scipy method, curve_fit, from optimise that invokes a non-linear sum of least squares method to calibrate our model to the given data. Using the standard deviations mentioned above as our sigma value for curve_fit and bounds of [0, 0, -∞] to [∞, ∞, ∞] for pressure and no bounds for temperature, we generated our best-fit parameters and their respective covariance matrix according to the local minimum associated with the starting values of [0.15 0.12 0.6] for pressure and [0.000001, 0.08] for temperature. Then using the covariance matrix as an input to np.random.multivariate_normal, we obtained a set of 100 adequality fitting models that show the range of pressure and temperature values.

![image](https://user-images.githubusercontent.com/85419997/174466479-5b36bab4-b5a1-486e-bd91-96632fc00867.png)
<br><i>Figure 15: Model forecast for pressure  for four scenarios. sampeled from a covariance matric generated by parameter association according by curve_fit</i>

![image](https://user-images.githubusercontent.com/85419997/174466489-8619cb05-78c2-4759-baa9-675aa816a157.png)
<br><i>Figure 16: Model forecast for temperature  for four scenarios. sampeled from a covariance matric generated by parameter association according by curve_fit</i>

In regards to pressure:
The 95% confidence interval is: 
-	Maintain [0.02643 	0.2832]
-	Stop [0.03958 		0.04266]
-	Double [0.01309 	0.0142]
-	Halve [0.03302 		0.03548]

In regards to temperature:
The 95% confidence interval is:
-	Maintain [145.08	145.12]
-	Stop [145.81		145.83]
-	Double [144.37		144.41]
-	Halve [145.44		145.48]

A more robust model would separate the cavities into smaller separate control volumes that interact with each other. While the overall reservoir pressure is above the threshold, local compartments could drop below the critical level and have a more significant effect on surface features 

Another shortcoming of our model is that extraction and injection rates are bundled into a single term. Separating these two rates into different terms would allow them to act independently of each other with varying parameter strengths, allowing for a better fit.

# 11.	Recommend?
The pressure and temperature change in the Rotorua Geothermal System have been modelled till 2050, subject to four different operation scenarios. Three of these scenarios (Halve, maintain, stop the current operation rate) maintained reservoir pressure above the 0.02 bar limit. Pressures below this limit have an adverse effect on surface features such as eruptions from the Waikite Geyser. The fourth scenario (double the current operation rate) drops the reservoir pressure to below 0 bar. 

In regards to temperature, three of these scenarios (Halve, maintain, stop the current operation rate) continued to increase in temperature while the fourth scenario (double the current operation rate) see the temperature start to slowly decline, which could have an adverse effect on surface features such as eruptions from the Waikite Geyser.

Suppose the intent is to avoid these impacts. In that case, reducing future operations to half production is recommended as this allows the Rotorua Geothermal System to recover faster while still allowing for the economic benefits of geothermal generation. 

The forecast above is specific to the assumptions we made and the scenarios we specified. Any deviation from this may show significant inconsistency with our predictions. For instance, if production was to be reduced by a third, there is no guarantee that the pressure would stay above the recommended limit. Furthermore, our model simplifies the reinjection of cold water and does not consider the interactions with surrounding rock stress states. In future studies, it would be better to a less investigate whether these had an impact on the pressure and temperature and hence the surface features.  


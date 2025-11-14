# Should-decommissioned-fossil-fuel-power-plants-be-preserved-as-industrial-heritage-
dataset and others
Data source
This dataset draws on a survey conducted within the framework of the INNOREJUST project. The full questionnaire and complete dataset can be consulted in:
García Docampo M, Sanz-Hernández A, Pérez-Sindín XS, Alonso-Domínguez À, Marcos Santiago M, Rodríguez Pacios A. Database and primary report related to the survey “The energy transition in Spain: Public opinion” [Dataset]. 2025. https://www.researchgate.net/publication/392324223
Compared with the original release, the present dataset includes only the variables relevant to our study. Each respondent can be linked back to the original dataset via the variable id.
Fieldwork
•	National sample (CATI): Conducted by the professional survey firm EDESGA (https://www.edesga.com/).
•	As Pontes (face-to-face): Conducted by EDESGA in collaboration with INZAREDE (https://www.inzarede.com/). Both firms used the same protocol and questionnaire software.
Replication
Analyses were run in R (RStudio). We share the R code used to reproduce the Principal Component Analysis (PCA) and the visualizations. The dataset is provided in CSV and SPSS (.sav) formats.
Notes on data quality
•	Missing values were imputed with subgroup weighted means.
•	When analyzing As Pontes, use WeightAsPontes; when analyzing Spain, use WeightSpain.
•	The As Pontes and national samples should be analyzed separately, even though they are stored in the same file.
Variables in the dataset
Below we provide each variable’s name, a short description, and—when applicable—the codes and response categories.
Survey – Survey sample
1 = Spain
2 = As Pontes
id – Respondent ID (numeric).
Municipality – Municipality of residence (string; 375 municipalities in data).
Code_ine – Official INE municipality code.
Province – Province of residence (54 provinces observed).
CCAA – Autonomous Community of residence
1 = Galicia
2 = Asturias
3 = Cantabria
4 = Basque Country
5 = Navarra
6 = La Rioja
7 = Aragón
8 = Madrid
9 = Castile and León
10 = Castile-La Mancha
11 = Extremadura
12 = Catalonia
13 = Valencian Community
14 = Balearic Islands
15 = Andalusia
16 = Region of Murcia
17 = Ceuta
18 = Melilla
19 = Canary Islands
Reindustrialization – Support for industrial redevelopment
0–10 = Scale values
11 = Don’t know / No answer
Renaturalization – Support for ecological restoration
0–10 = Scale values
11 = Don’t know / No answer
Preservation – Support for heritage reuse
0–10 = Scale values
11 = Don’t know / No answer
WeightAsPontes – Survey weight for As Pontes sample (continuous).
WeightSpain – Survey weight for Spain sample (continuous).
NearestPowerPlantName – Name of nearest coal-fired power plant
•	Aboño
•	Andorra
•	As Pontes
•	Carboneras
•	Cercs
•	Escatrón
•	Escucha
•	La Robla
•	Los Barrios (note: actually in Cádiz, mis-coded under Tenerife in dataset)
•	Puertollano
•	Velilla
•	Compostilla II
•	La Pereda
•	Meirama
•	Narcea
•	Puente Nuevo
•	Soto de Ribera
PowerPlantCluster – province in which the nearest power plant is located
•	Almería
•	Asturias
•	Ciudad Real
•	A Coruña
•	Córdoba
•	León
•	Lleida
•	Palencia
•	Cádiz
•	Teruel
Dist_to_nearest_plant – Distance to nearest power plant in kilometers (range ≈ 1.38 to 1416.73).
<img width="432" height="643" alt="image" src="https://github.com/user-attachments/assets/a6429197-f0ca-4a4c-8c7c-e7dc6ae26fac" />

#Denne Bachelor indeholder følgende:

1. Replikation_af_artiklen.md: Indeholder koden til replikationen af resultaterne fra Acerbi og Szekelys artikel (2014).

2. S&P_500_Historical_Data.csv: Indeholder de empiriske data, der anvendes i projektet. Datasættet dækker perioden 1/1-2020 – 31/1-2025.

3. Implementering_af_data.md: Indeholder koden til implementeringen af data fra S&P 500-indekset. Vi opdeler datasættet i trænings- og testdata og fitter en AR(1)-GARCH(1,1)-model til træningsdataene.

4. Sensitivitetsanalyse.md: Indeholder koden til vores sensitivitetsanalyse, hvor vi undersøger robustheden af fire teststørrelser: General CC udviklet af Nolde og Ziegel, Intercept ESR foreslået af Bayer og Dimitriadis samt teststørrelserne Z1 og Z2 fra Acerbi og Szekely.
   
5. Poweranalyse.md: Tester de fire ovenstående teststørrelsers styrke (power) under både historisk simulering og filtreret historisk simulering.

6. Test_på_empirisk_data.md: Indeholder koden til test af vores model på empiriske data fra S&P 500-indekset (testdata). Datasættet dækker perioden 1/2-2024 – 31/1-2025. Vi analyserer, hvor godt vores AR(1)-GARCH(1,1) model passer til empiriske data.

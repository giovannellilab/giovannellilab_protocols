![wetbench](https://img.shields.io/badge/TYPE-wet_bench-brigthgreen)
![author](https://img.shields.io/badge/AUTHOR-Donato_Giovannelli-ad7fa8)
![created](https://img.shields.io/badge/created-25092020-lightgray)
[![giovannellilab](https://img.shields.io/badge/BY-Giovannelli_Lab-blue)](http://dgiovannelli.github.io)

 <img src="https://dgiovannelli.github.io//images/logopic/giovannellilab.png" width="100 px">

# Determination of Spore abundance in soils/sediments using DPA concentrations as a proxy

**OBJ**: Estimate the amount of endospores in a sediment/soil sample using Dipicolonic acid as a proxy for endospore abundance

>***Duration: not yet assessed***

> ***Adapted from: Fitchel et al. 2007 J Microbiol Methods***

#### Materials
- 15 ml falcon tubes (autoclave safe)
- 2.2 ml tubes
- Ice bath
- 0.22 µm cellulose Anotop filter
- P5000, P1000, P200 tips
- 60°C Oven

#### Equipment
- P1000
- P5000
- P200
- Centrifuge for 15 ml falcon tubes
- Autoclave
- Spectrofluorimeter
- Quartz cuvette


#### Solutions
- KCl 50 mM
- Sodium bisulphate buffer 50 mmol/l, pH 1.2
- Tb3Cl in 400 mM sodium acetate buffer pH 5

## Procedure
1. Number and pre-weight (record the weight on the DPA spreadsheet) enough 15 ml falcon tubes + 2 blanks
2. Weight 2 g of wet sediments into the falcon
3. Add 2 mL KCl (50 mM) and vortex
4. Centrifuge 5 min @ X 7,000 g and discard supernatant. This removes readily exchangeable calcium
5. Add a second 2 mL KCl (50 mM) aliquot and vortex
6. Centrifuge 5 min @ X 7,000 g and discard supernatant
7. Add 3 ml sodium bisulphate buffer (50 mmol/l, pH 1.2) (Sodium acetate 500 mM pH 5 also possible)
8. Autoclave 30 min @ 121°C
9. Cool on ice
10. Centrifuge 5 min @ X 5,000 g
11. Filter the supernatant with a 0.22 µm cellulose acetate filter into a new 2.2 ml vial
12. Dry the sediment pellet for 24 h @ 60°C and weight the dry pellet. The next day record the pellet weight
13. Set the fluorescence blank using 2 ml of sodium bisulphate buffer (50 mmol/l, pH 1.2)
14. Load 2 ml of the supernatant into the quartz cuvette (read the two reagent blank first!)
15. Read on spectrofluorimeter (@ ex. 270 nm and em. 548 nm) and record the background fluorescence
16. Add XXX µl of Tb3Cl (final 30 µM in 400 mM sodium acetate buffer pH 5) to the cuvette
17. Briefly mix by pipetting or inverting several time the cuvette
18. Read again (@ ex. 274 nm and em. 548 nm) and record the terbium-DPA complex fluorescence
19. Discard and wash the cuvette with ddH2O water
20. Repeat until all sample have been read for background and terbium-DPA fluorescence
21. Calculate the equivalent endospore abundance using the formula below and the appropriate spreadsheet

## Expected results (quantitative information, graphics, images)
Expected results are 2 spectrofluirmetric reads per sample (background and terbium-DPA complex) plus two reagent blanks reads. There can be converted, together with the sediment dry weight of each sample into the DPA concentration using the calibration curve and into endospore abundance using a conversion factor as outlined below.

## Common problems, troubleshooting and solutions
- None to date

## Calculations
Calculations take into account the reagent blank fluorescence, the background fluorescence, the dilution factor, a recovery correction factor and the starting sediment dry weight. All this is converted to DPA abundance within a range of X-XXX µM of DPA.
A correction factor of 1.2 was used on the obtained DPA concentration to account for known recovery efficiency from sediments (Fitchel et al., 2007 J Microbiol Methods). A DPA to spore conversion factor of 2.24 x 10^-16 mol/endospore (Fitchel et al. 2007 FEMS Microbiol Ecol) was used, although the original estimate suggest a variation of plus/minus 0.69x 10^-16 mol. Additionally, an increased concentration in DPA per endospore was observed fro spores sporulating at 40°C vs 25°C with a DPA per endospore content doubling at higher sporulation temperatures. Although we have no evidence of this effect in natural samples, nor it is possible to know the real temperature of sporulation, this should be taken into account in future work by applying a temperature dependent correction factor.

$$
DPA=(((m\times((Fs-Fb)-Fk)\timesdil)+q)\timesRc)\timesf
Ea=\frac{DPA\timesEc}{Sw}
$$

where:
DPA = DPA concentration in µmol/l
Ea = Endospores abundance in endospores/g dry sediment
Fk = Fluorescence blank reagents
Fb = Background fluorescence before Tb3Cl addition
Fs = Sample fluorescence after Tb3Cl addition
Sw = Sediment dry weight
dil = dilution factor (if used)
f = Fractional factor correcting for protocol volume handling
Rc = Recovery correction
m = slope of the calibration curve
q = intercept of the calibration curve
Ec = Endospore DPA content

A Open Document Spreadsheet is attached to this protocol to perform the calculations.

> Written with [StackEdit](https://stackedit.io/)

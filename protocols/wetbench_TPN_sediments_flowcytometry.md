![wetbench](https://img.shields.io/badge/TYPE-wet_bench-brigthgreen)
![author](https://img.shields.io/badge/AUTHOR-Donato_Giovannelli-ad7fa8)
![created](https://img.shields.io/badge/created-25092020-lightgray)
[![giovannellilab](https://img.shields.io/badge/BY-Giovannelli_Lab-blue)](http://dgiovannelli.github.io)

 <img src="https://dgiovannelli.github.io//images/logopic/giovannellilab.png" width="100 px">

# Determination of cell abundance in soils/sediments using flow cytometry

**OBJ**: Determination of total prokaryotic numbers (TPN) and biomass (TPB) in formaldehyde-fixed sediment samples using flow cytometry. Two different variation are proposed, one with direct staining, the other filter based. Please assess which of the two version is more suitable for your sample type. Asses also if the sample after treatment appears to be too dirty for direct flow cytometry counting.

>***Duration: not yet assessed***

> ***Adapted from: Morono et al., 2013 Env Microbiol. See also Amalfitano protocols***

#### Materials
- list all the material

#### Equipment
- list all the qeuipment

#### Solutions
- list the necessary solutions. Be sure to include the concentration, pH and notes on the preparation or and any reference to the protocols to prepare them

## Procedure
1. Mix vigorously the stored sample to make an homogeneous slurry
2. Make sure that the slurry is homogeneous. This is a critical step!
3. Transfer 250 µl of slurry (using a cut 1 ml tip or a small spoon) to a new 2 ml tubes
4. Add 200 µl of detachment solution (100 mM EDTA, 100 mM Na-pyrophosphate, 1% v/v Tween 80)
5. Add 50 µl methanol
6. Shake on vortex (with vortex adapter) 60' @ 500 rpm
7. Pulse sonicate in water bath 1' @ 20 W
8. Transfer 50 µl of cell-detached sediments with 500 µl of TE buffer in a 40 µm cell mesh strainer
9. Sieve to remove large mineral particles

#### Direct staining version of the protocol (direct staining version)
10. Place the sieved cell-detached sediment suspension in a flow cytometry tube
11. Add 1 ml of TE buffer
12. Add 100 µl of Sybr Green I and stain for 5'
13. Count using the Accuri C6 using both the green (525 nm) and red (695 nm) fluorescence. Count 10,000 events within the high red fluorescence gate. Cells have a higher green fluorescence compared to the SYBR-SPAM (extracellular organic carbon)

#### Filter-based staining version of the protocol (filter staining version)
10. Place the sieved cell-detached sediment suspension onto a 0.22 µm Anodisc filter
11. Wash with 1 ml TE buffer
12. Stain the filter with 100 µl of Sybr Green I and stain for 5'
13. Place the filter in a 15 ml falcon tube with 2 ml TE buffer
14. Sonicate 20" @ 20 W
15. Transfer the supernatant to a clean 15 ml falcon tube
16. Add 1 ml of TE buffer to the filter and sonicate again
17. Combine the supernatants
18. Centrifuge 10' @ 6,000g
19. Remove 2 ml of supernatant to reduce the volume
20. Resuspend cells in the remaining volume and transfer to a flow cytometry tube
21. Count using the Accuri C6 using both the green (525 nm) and red (695 nm) fluorescence. Count 10,000 events within the high red fluorescence gate. Cells have a higher green fluorescence compared to the SYBR-SPAM (extracellular organic carbon)

## Expected results (quantitative information, graphics, images)
Expected results are the flow cytometry event count exported as CSV file, the raw machine data, the number of counted events within the gated area, the total volume counted, the sediment dry weight. TPN abundances can be calculated with the formula below and the associated spreadsheet file.

## Common problems, troubleshooting and solutions
- Low red fluorescence of cells. this is due to preferential absorption of the SYBR dye to the sediment organic matrix. Increase the SYBR concentrations either restaining the sample a second time or using a fresh stock of SYBR at a higher working concentration.

## Calculations
The total prokaryotic abundance can be calculated as follow:

$$
TPN=\frac{\frac{Ev}{Vp}\times{dil}\times{10}}{Sw}
TPB=(TPN\times{Cc})\times{10^9}
$$

where:
TPN = Total Prokaryotic Numbers in n. cell / g of dry sediment
TPB = Total Prokaryotic Biomass in µg C / g of dry sediment
E_v = Number of events in the gate
V_p = Volume processed in the cytometer
V_t = Total extraction volume
dil = dilution factor
f = Fractional factor correcting for protocol volume handling
S_w = Sediment dry weight
C_c = Cell carbon content in fg C / cell





To obtain replicates, the sytometry run can be divided in time to obtain subreplicae.

> Written with [[StackEdit](https://stackedit.io/)

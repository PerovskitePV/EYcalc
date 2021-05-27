# Energy Yield Calculator

<p align="center"><img src="https://raw.githubusercontent.com/wiki/PerovskitePV/EYcalc/Logo.png"></p>

This software aims to simulate the energy yield of single-junction and multi-junction solar cells. In contrast to the power conversion efficiency (PCE), the energy yield (EY) accounts for environmental conditions, such as constantly changing irradiation conditions or the ambient temperature.

It allows rapid simulation of complex architectures and was developed with the aim to handle textured perovskite-based multi-junction devices. However, it is possible to simulate any combination of thin-film architecture with incoherent photovoltaic materials (e.g., crystalline silicon).

By making use of pre-simulated textures (e.g., inverted pyramids, regular upright pyramids, random pyramids) by geometrical ray tracing, any incoherent interface within the architecture can also be textured. 

The software is available as source code and a simple to use graphical user interface (GUI), which requires either a MATLAB (>R2017a) installation or the MATLAB runtime.

### Basic Features

The basic features of the **EYcalc** are:

* Spectral and angular-resolved realistic irradiance data (from 1020 locations in the USA) is used
* A simple cloud model is used to adjust the diffuse irradiation
* Fast optical simulations, by combining the transfer matrix method and geometric ray tracing
* Optics can handle arbitrary combinations of thin (coherent) and thick (incoherent) layers, which also can be textured
* Single- and multi-junction solar cells can be simulated
* No limitation on the number of absorbers
* Energy yield is computed for different electrical interconnection schemes (e.g. 2T, 3T, 4T)
* Energy yield can be derived for constant tilt (and constant rotation) angle
* Energy yield can be derived for various tracking algorithms (e.g. 1-axis, 2-axis)
* Bifacial solar cells can be simulated
* Albedo can be considered by choosing one out of 3400 spectra of natural and man-made materials from the [ECOSTRESS spectral library](https://speclib.jpl.nasa.gov/)

### Modular framework

The software is divided into individual modules, which handle the irradiation, optics, electrics and energy yield simulations. Those modules can also be operated independently (e.g., calculate the reflectance, transmittance, absorptance of a solar cell architecture).

[![](https://mermaid.ink/img/eyJjb2RlIjoiZ3JhcGggTFIgXG5cdFRNWTNbTWV0ZW9yb2xvZ2ljYWwgRGF0YSBUTVkzXSAtLT4gSXJyW0lycmFkaWFuY2UgTW9kdWxlXSAtLT4gRVlDW0VuZXJneSBZaWVsZCBDb3JlIE1vZHVsZV0gLS0-IEVZW0VuZXJneSBZaWVsZF1cblx0RFtEZXZpY2UgQXJjaGl0ZWN0dXJlXSAtLT4gT1tPcHRpY3MgTW9kdWxlXSAtLT4gRVlDXG5cdEVZQyAtLT4gRVtFbGVjdHJpY3MgTW9kdWxlXVxuXHRFIC0tPiBFWUNcblx0RVBbRWxlY3RyaWNhbCBQcm9wZXJ0aWVzXSAtLT4gRSIsIm1lcm1haWQiOnsidGhlbWUiOiJkZWZhdWx0In0sInVwZGF0ZUVkaXRvciI6ZmFsc2V9)](https://mermaid-js.github.io/mermaid-live-editor/#/edit/eyJjb2RlIjoiZ3JhcGggTFIgXG5cdFRNWTNbTWV0ZW9yb2xvZ2ljYWwgRGF0YSBUTVkzXSAtLT4gSXJyW0lycmFkaWFuY2UgTW9kdWxlXSAtLT4gRVlDW0VuZXJneSBZaWVsZCBDb3JlIE1vZHVsZV0gLS0-IEVZW0VuZXJneSBZaWVsZF1cblx0RFtEZXZpY2UgQXJjaGl0ZWN0dXJlXSAtLT4gT1tPcHRpY3MgTW9kdWxlXSAtLT4gRVlDXG5cdEVZQyAtLT4gRVtFbGVjdHJpY3MgTW9kdWxlXVxuXHRFIC0tPiBFWUNcblx0RVBbRWxlY3RyaWNhbCBQcm9wZXJ0aWVzXSAtLT4gRSIsIm1lcm1haWQiOnsidGhlbWUiOiJkZWZhdWx0In0sInVwZGF0ZUVkaXRvciI6ZmFsc2V9)

* The [**Irradiation Module**](https://github.com/PerovskitePV/EYcalc/wiki/Irradiance-Module) calculates the spectral and angular-resolved irradiance over the course of one year with a temporal resolution of one hour by applying [SMARTS](https://www.nrel.gov/grid/solar-resource/smarts-register.html) to typical meteorological year ([TMY3](https://nsrdb.nrel.gov/data-sets/archives.html)) data of locations in various climatic zones. A simple model is employed to account for cloud coverage such that realistic direct and diffuse irradiance are derived.
* The [**Optics Module**](https://github.com/PerovskitePV/EYcalc/wiki/Optics-Module) rapidly calculates the spectral and angular-resolved absorptance of the non-simplified architecture of multi-junction solar cells. It is able to handle multiple planar and textured interfaces with coherent and incoherent light propagation by combining transfer matrix method (TMM) and geometrical ray-tracing.
* The [**Electrical Module**](https://github.com/PerovskitePV/EYcalc/wiki/Electrics-Module) determines the temperature-dependent current density-voltage (*J*-*V*) characteristics accounting for series and shunt resistances for a given short-circuit current density (*J*<sub>SC</sub>) of the sub-cells forming the multi-junction in either a 2T-, 3T- or 4T-configuration. Furthermore, the maximum power point is determined to calculate the power output of the multi-junction solar module.
* The [**Energy Yield Core Module**](https://github.com/PerovskitePV/EYcalc/wiki/Energy-Yield-Module) calculates the EY over the course of one year of the sub-cells depending on their orientation (rotation and/or tilt of the module) and location. The EY is computed by combining the spectral and angular resolved solar irradiation (with or without albedo), the absorptance of the multi-junction solar cell and the electrical properties.

### Credits

This software project was initiated by **[Ulrich W. Paetzold](mailto:ulrich.paetzold@kit.edu?subject=[GitHub]%20Question%20on%20Energy%20Yield%20Software)**. The **code development** was driven by:

* **Jonathan Lehr** (electrics module, albedo)
* **Malte Langenhorst** (optics module, irradiance module)
* **Raphael Schmager** (energy yield core, irradiance module, optics module, electrics module, GUI)
* **Fabrizio Gota** (numerical modelling on 3T interconnection, optics module)

The financial support by the following **projects and grants** is gratefully acknowledged:

- [PERCISTAND](https://percistand.eu/en) (funding code: 850937), European Union's Horizon 2020 research and innovation programme
- Helmholtz Young Investigator Group of U. Paetzold (funding code: VH-NG-1148), [Helmholtz Association](https://www.helmholtz.de/)
- [PEROSEED](https://www.helmholtz-berlin.de/projects/peroseed/index_en.html) (funding code: ZT-0024), [Helmholtz Association](https://www.helmholtz.de/)
- CAPITANO (funding code: 03EE1038B), [Federal Ministry for Economic Affairs and Energy](https://www.bmwi.de/)
- 27Plus6 (funding code: 03EE1056B), [Federal Ministry for Economic Affairs and Energy](https://www.bmwi.de/)

This software uses codes and data from **other programmers and resources**:

* Parts of the transfer matrix code is taken from [Steven Byrnes](https://github.com/sbyrnes321/)
* Matlab implementation of the [NREL solar position algorithm](https://doi.org/10.1016/j.solener.2003.12.003) by [Vincent Roy](https://de.mathworks.com/matlabcentral/fileexchange/5430-sun-azimuth-data) 
* Logarithmic Lambert W function from [Michael](https://www.mathworks.com/matlabcentral/fileexchange/57239-lambert-w-function-logarithmic-input)
* The [SMARTS](https://www.nrel.gov/grid/solar-resource/smarts-register.html) from Dr. Christian A. Gueymard [see also](https://www.solarconsultingservices.com/smarts.php)
* The [TMY3](https://nsrdb.nrel.gov/data-sets/archives.html) data from the National Solar Radiation Database
* [Reference Air Mass 1.5 Spectra](https://www.nrel.gov/grid/solar-resource/spectra-am1.5.html)
* [ECOSTRESS spectral library](https://speclib.jpl.nasa.gov/) for albedo

### Getting started

Download and extract the project. Open the `main.m`, which contains all definitions and settings to calculate the energy yield (EY) of a perovskite/c-Si multi-junction solar cell. In order to start your simulations, you need to get and add some external files, like the SAMRTS code or the TMY3 data. Please check out our [wiki page](https://github.com/PerovskitePV/EYcalc/wiki) for some help. You'll find a detailed description of each of the modules as well as a [guide](https://github.com/PerovskitePV/EYcalc/wiki/Setup) for setting up the required external files.   

### Contributing

If you want to contribute to this project and make it better, your help is very welcome! 

### Contact

For any questions regarding the software, please contact [Ulrich W. Paetzold](mailto:ulrich.paetzold@kit.edu?subject=[GitHub]%20Question%20on%20EYcalc). See also [here](https://www.lti.kit.edu/mitarbeiter_7254.php).

### Citing

If you use our software or parts of it in the current or a modified version, you are obliged to provide proper attribution. This can be to our paper describing the software:

* R. Schmager and M. Langenhorst et al., Methodology of energy yield modelling of perovskite-based multi-junction photovoltaics, Opt. Express. (2019). [doi:10.1364/oe.27.00a507](https://doi.org/10.1364/OE.27.00A507).

or to this code directly:

* EYcalc - Energy yield calculator for multi-junction solar modules with realistic irradiance data and textured interfaces. (2021). [doi.org/10.5281/zenodo.4696257](https://doi.org/10.5281/zenodo.4696257).

### License

This software is licensed under the [GPLv3](https://www.gnu.org/licenses/gpl-3.0.html) license. © 2021 EYcalc -
Ulrich W. Paetzold, Raphael Schmager, Malte Langenhorst, Jonathan Lehr, Fabrizio Gota

Interested in a sublicense agreement to use EYcalc in a non-free/restrictive environment? Contact [Ulrich W. Paetzold](mailto:ulrich.paetzold@kit.edu?subject=[GitHub]%20Question%20on%20EYcalc)!

### Further reading

This energy yield software has been used in the following publications:

* M. De Bastiani et al., Efficient bifacial monolithic perovskite/silicon tandem solar cells via bandgap engineering, Nature Energy. (2021). [doi.org/10.1038/s41560-020-00756-8](https://doi.org/10.1038/s41560-020-00756-8).

* J. Lehr et al., Numerical study on the angular light trapping of the energy yield of organic solar cells with an optical cavity, Opt. Express. (2020). [doi.org/10.1364/OE.404969](https://doi.org/10.1364/OE.404969).

* F. Gota et al., Energy Yield Advantages of Three-Terminal Perovskite-Silicon Tandem Photovoltaics, Joule, (2020). [doi.org/10.1016/j.joule.2020.08.021](https://doi.org/10.1016/j.joule.2020.08.021).

* J. Lehr et al., Energy yield of bifacial textured perovskite/silicon tandem photovoltaic modules, Sol.
  Energy Mater. Sol. Cells. (2020). [doi:10.1016/j.solmat.2019.110367](https://doi.org/10.1016/j.solmat.2019.110367).

* R. Schmager et al., Methodology of energy yield modelling of perovskite-based multi-junction
  photovoltaics, Opt. Express. (2019). [doi:10.1364/oe.27.00a507](https://doi.org/10.1364/OE.27.00A507).

* M. Langenhorst et al., Energy yield of all thin-film perovskite/CIGS tandem solar modules, Prog.
  Photovoltaics Res. Appl. (2019). [doi:10.1002/pip.3091](https://doi.org/10.1002/pip.3091). 

* J. Lehr et al., Energy yield modelling of perovskite/silicon two-terminal tandem PV modules with flat
  and textured interfaces, Sustain. Energy Fuels. (2018). [doi:10.1039/c8se00465j](https://doi.org/10.1039/C8SE00465J).
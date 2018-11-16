## GEARS
by Jake Alan Pitt and Julio R. Banga

<a href="https://doi.org/10.5281/zenodo.1420464"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.1420464.svg" alt="DOI"></a>

GEARS (Global parameter Estimation with Automated Regularisation via Sampling) is a Matlab toolbox for parameter estimation in nonlinear dynamic models composed of deterministic ordinary differential equations (ODEs). GEARS is based on the combination of three main strategies: (i) global optimisation (to avoid convergence to local solutions), (ii) reduction of the search space (i.e. tighter bounds on parameters), (iii) regularised estimation, a strategy used to handle overfitting (i.e. fitting the noise rather than the signal). As a result, GEARS can allow the user to avoid both underfitting and overfitting problems while requiring minimum supervision from the user. These capabilities are especially useful when calibrating ODE-based models with highly nonlinear and flexible dynamics (e.g. models of biological oscillators).

## Features

<ul>
<li>Unsupervised parameter estimation </li>
<ul>
<li>Efficient Global optimisation</li>
<li>Automated regularisation tuning</li>
<li>Parameter bounding</li>
</ul>
<li>Post-fit analyses</li>
<ul>
<li>Normalised root mean square error (NRMSE)</li>
<li>R<sup>2</sup> test</li>
<li>chi<sup>2</sup> test</li>
<li>Parameter uncertainty (using Fisher information)</li>
<li>Correlation matrix (using Fisher information)</li>
<li>Active bounds</li>
</ul>
<li>Post-fit plotting</li>
<ul>
<li>Trajectories for Fitting and cross-validation</li>
<li>Parameter space samples</li>
<li>Visualisation of parameter bound reduction</li>
<li>Residuals</li>
<li>Predictions vs measurements</li>
<li>Trajectory uncertainty (via Fisher information)</li>
<li>Convergence curves</li>
</ul>
<li>Results reports</li>
<ul>
<li>html markup of results</li>
<li>xls markup of results</li>
<li>Combined pdf report of all plotted Figures</li>
</ul>
</ul>

## Requirements

<p> GEARS runs on Matlab R2015b or later and is multi-platform (Windows and Linux). Both the optimisation and symbolic mathematics Matlab toolboxes are required to run GEARS. </p>

<p>GEARS requires that the <a href="http://icb-dcm.github.io/AMICI/">AMICI package</a> has been correctly installed.</p>

<p>Optionally, users can use <a href="https://www.ghostscript.com">Ghostscript</a> for the exportation of figure reports.</p>

## Citations 

<p>If you use this toolbox and publish the results, please cite it with the following references.</p>

<p> Regarding the methodology: <br>
Pitt, J.A. and Banga, J.R. (2018) Parameter estimation in models of biological oscillators:
an automated regularised estimation approach. Submitted. </p>

<p> Regarding the software:<br>
Pitt, J.A. and Banga, J.R. (2018) GEARS - a toolbox for Global parameter Estimation
with Automated Regularisation via Sampling. <br>
doi: 10.5281/zenodo.1420465</p>

## Documentation

GEARS' documentation is available in PDF format. Said documentation gives a detailed description of how to install and use the GEARS package with illustrative examples.

## License 

GEARS is distributed under the  <a href="http://www.gnu.org/licenses/gpl.html">GNU General Public License version 3 (GPL v3)</a>. Copyright 2018 Jake Alan Pitt and Julio R. Banga.

## Acknowledgements 

GEARS was developed at (Bio)Process Engineering group - IIM-CSIC, Spanish National Research Council c/Eduardo Cabello, 6. 36208, Vigo (Spain). This research received funding from the European Union's Horizon 2020 research and
innovation program under grant agreement No 675585 (MSCA ITN \SyMBioSys") and
from the Spanish MINECO/FEDER project SYNBIOCONTROL (DPI2017-82896-C2-2-R). Jake Alan Pitt is a Marie Sklodowska-Curie Early Stage Researcher under the supervision of Prof. Julio R. Banga.

## Support

<p> Please check the document and examples carefully before contacting the authors. Please
pay special attention to section 9 of the documentation. If needed, they can be contacted by
email at: </p>
Jake Alan Pitt - jp00191.su@gmail.com <br>
Julio R. Banga - julio@iim.csic.es.

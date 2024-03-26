# CISER GerryChain Workshop

Home for the materials (notebooks, data, and commentary) associated to the Spring 2024 WSU CISER Workshop on Computational Redistricting. This tutorial is intended to provide a starting point for exploring computational redistricting in Python and doesn't assume any specific programming background. The intended outline is approximately: 

* How do Jupyter notebooks work in Python?
* Data Processing for Redistricting
* Metrics on Redistricting Plans
* The MCMC Ensemble Method

with the goal of providing lots of examples and opportunities for attendees to experiment.  

## Installation Instructions
The Python Install .pdf above provides some links and instructions for getting started with Python on your own computer. If you would prefer not to install anything now, you can still follow along by either using the directions in that same .pdf for opening the notebooks with one of the cloud services or by simpling clicking the corresponding links above. 

## Notebooks
The main material we will be discussing in the tutorial is contained in the four Jupyter notebooks presented above. 

*  **0\_Jupyter\_Basics.ipynb** This notebook walks through the basic usage of Jupyter notebooks and how to interact with Python cells. It also includes examples of  basic arithmetic, variable assignment, data types, data structures,  and imports work in Python as well as a couple of exercises to check your understanding. If you already have some familiarity with Python, this one is safe to skip :)
* **1\_Dual\_Graphs\_and\_MAUP.ipynb** This notebook explores the basic properties of dual graphs (using the networkx package) for representing geospatial data (using geopandas) for redistricting.  We will also use the MAUP package to move demographic and partisan data between different levels of resolution and assign census units to plans. 
* **2\_Compactness.ipynb** In this notebook we will examine continuous and discrete geographic measures on districting plans. 
* **3\_Partisan\_Symmetry.ipynb** In this notebook we use the Partition object from the GerryChain package to analyze the partisan symmetry of redistricting plans. 
* **4\_Ensembles.ipynb** In this final notebook we will put everything together, using tree-based MCMC methods to generate ensembles of redistricting plans. 


<table>
  <tr><td>
<img src="https://github.com/drdeford/CISER_GerryChain/blob/main/Animations/LWAR.gif" width=300>
    </td><td>
<img src="https://github.com/drdeford/CISER_GerryChain/blob/main/Animations/WA4.gif" width=300>
        </td><td>
<img src="https://github.com/drdeford/CISER_GerryChain/blob/main/Animations/forest_lifted_walk.gif" width=300>
    </td><td>
<img src="https://github.com/drdeford/CISER_GerryChain/blob/main/Animations/NH_Two_Edge_Steps.gif" width=300>
    </td></tr>
    <tr><td>
<img src="https://github.com/drdeford/CISER_GerryChain/blob/main/Animations/lifted_gif.gif" width=300>
    </td><td>
<img src="https://github.com/drdeford/CISER_GerryChain/blob/main/Animations/new_wa_gif.gif" width=300>
    </td><td>
<img src="https://github.com/drdeford/CISER_GerryChain/blob/main/Animations/space10.gif" width=300>

</td><td>
<img src="https://github.com/drdeford/CISER_GerryChain/blob/main/Animations/ID_plans_sld.gif" width=300>

</td></tr>
</table>

## Other Resources

A couple of hours isn't nearly enough time to explore all of the facets of computational redistricting, so here are some good starting places if you are interested in learning more: 
* If you are interested in trying out the GerryChain package for analyzing districting plans, you can find:
  * General information and background: htps://tinyurl.com/gerrytalk
  * Mathematical background (some typos!): http://math.wsu.edu/faculty/ddeford/mcmc_intro.php
  * Templates and notebooks here: https://github.com/drdeford/GerryChain-Templates 
  * A bootcamp-style introduction here: https://github.com/vrdi/GerryChain-BootCamp
* The NetworkX documentation is always a good starting spot and contains numerous examples in addition to the descriptions of functions and parameters: https://networkx.org/documentation/stable/index.html
* I taught a graduate class on Computational Tools for Complex Networks at WSU last semester and the course repository (with additional notebooks and exposition documents) can be found here: https://github.com/drdeford/Math_581_05 The notebooks don't have as much expository material as the ones for this tutorial session but do include more detailed computations.
* As a part of the 2019 Voting Rights Data Institute, I taught a week of breakout sessions about networks with materials here: https://github.com/vrdi/Networks-Breakout


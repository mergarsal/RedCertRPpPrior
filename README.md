# Certifiable solvers for the relative pose problem with known gravity vector

Estimate the relative pose between two calibrated cameras 
with known gravity vector. 
DIfferent formulations of the problem 
are used to certify the solution
Refer to our paper [HERE](https://www.sciencedirect.com/science/article/pii/S0004370223000085?via%3Dihub) 
for more information.



**Authors:** 
[Mercedes Garcia-Salguero](https://mapir.isa.uma.es/mapirwebsite/?p=1718), 
[Javier Gonzalez-Jimenez](https://mapir.isa.uma.es/mapirwebsite/?p=1536)


**License:** [GNUv3](https://github.com/mergarsal/RedCertRPpPrior/blob/main/LICENSE)


If you use this code for your research, please cite:

```
@ARTICLE{,
    author = {Garcia-Salguero, Mercedes and Gonzalez-Jimenez, Javier},
     title = {Fast certifiable relative pose estimation with gravity prior},
   journal = {Artificial Intelligence},
      year = {2023},
       url = {http://mapir.isa.uma.es/papersrepo/2023/2023_mercedes_AI_priorRPp_doc.pdf},
       doi = {https://doi.org/10.1016/j.artint.2023.103862}
}
```



# Dependencies

The certifier requires : 
1. *Optimization* by D. Rosen. 
We use our fork

```
https://github.com/mergarsal/Optimization
```
2. *Iterative certifier* by us. 
```
https://github.com/mergarsal/QCQPIterCertifier
```

Use 
```
git submodule update --init --recursive
```
To download the dependency


## Build
```
git clone https://github.com/mergarsal/RedCertRPpPrior.git
cd RedCertRPpPrior

mkdir build & cd build 

cmake .. 

make -jX

```

The compiled examples should be inside the `bin` directory. Run: 
```
        ./bin/example_basic
```
 


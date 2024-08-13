# codeCollection
My source code collection, especially CFD codes

## From book
### FDM·FVM
* [エクセルでできる熱流体のシミュレーション](https://www.maruzen-publishing.co.jp/item/b304713.html), Excel
* [Fluid Flow Phenomena: A Numerical Toolkit](https://www.amazon.co.jp/Fluid-Flow-Phenomena-Applications-2006-01-01/dp/B01JXQ96NS), Fortran
* [乱流の数値シミュレーション](https://www.yokendo.com/books/9784842505268/), Fortran
* [流体解析の基礎](https://www.asakura.co.jp/detail.php?book_code=13111), Fortran
* [数値計算による流体力学-ポテンシャル流，層流，そして乱流へ-](https://www.coronasha.co.jp/np/isbn/9784339046519/), Fortran
* [Pythonで学ぶ流体力学の数値計算法](https://www.ohmsha.co.jp/book/9784274224706/), Python
* [数値流体解析の基礎 Visual C++とgnuplotによる圧縮性・非圧縮性流体解析](https://www.coronasha.co.jp/np/isbn/9784339046649/), C++
* [数値流体工学](https://www.utp.or.jp/book/b302011.html), Fortran
* [I do like CFD](http://ossanworld.com/cfdbooks/cfdbooks.html), Fortran
* [コンピュータによる流体力学](https://www.maruzen-publishing.co.jp/item/b294565.html), Fortran
* [計算流体力学-CIPマルチモーメント法による手法 CD-ROM付-](https://www.coronasha.co.jp/np/isbn/9784339045970/), C++
### FEM
* [詳解流れの数値計算-有限要素法による非圧縮性流体解析の基礎-](https://www.coronasha.co.jp/np/isbn/9784339046762/), Fortran
* [有限要素法による流れのシミュレーション](https://www.maruzen-publishing.co.jp/item/b295175.html), Fortran
### LBM
* [格子ボルツマン法](https://www.morikita.co.jp/books/mid/067661), C
* [格子ボルツマン法入門: 複雑境界および移動境界流れの数値計算法](https://www.maruzen-publishing.co.jp/item/b303613.html), C++
* [格子ボルツマン法・差分格子ボルツマン法](https://www.coronasha.co.jp/np/isbn/9784339046588/), Fortran
### MPS
* [粒子法入門](https://www.maruzen-publishing.co.jp/item/b294791.html), C++

## From web and repository
### FDM·FVM
* [SLFCFD](https://slfcfd.sourceforge.net/), C, 2-D and 3-D compressible flow and for 2-D incompressible flow.
* [nast2d](https://ins.uni-bonn.de/content/software-nast2d), C, Fortran, 2D, incompressible Navier-Stokes equations.
* [CADMAS-SURF](https://www.cdit.or.jp/program/cadmas.html), Fortran, wave/free surface flow.
* [PARIS](http://www.ida.upmc.fr/~zaleski/paris/), Fortran, multiphase flow solver with VOF and fronttracking methods.
* [surfer](http://www.ida.upmc.fr/~zaleski/codes/legacy_codes.html), Fortran, droplet and bubble legacy codes.
* [afid](https://stevensrjam.github.io/Website/afid.html), Fortran, canonical wall bounded turbulent flows.
* [pencil-code](http://pencil-code.nordita.org/), Fortran, compressible hydrodynamic flows with magnetic fields.
* [CFL3D](https://nasa.github.io/CFL3D/), Fortran, structured-grid, cell-centered, upwind-biased, Reynolds-averaged Navier-Stokes (RANS) code.
* [xcompact3d/Incompact3d](https://www.incompact3d.com/), Fortran, DNS/LES for turbulent flows.
* [CaNS](https://github.com/p-costa/CaNS), Fortran, canonical Navier-Stokes solver with FFT.
* [REEF3D](https://reef3d.wordpress.com/), C++, coastal, marine and hydraulic engineering flows.
* [gerris](http://gerris.dalembert.upmc.fr/), C, free software program for the solution of the partial differential equations describing fluid flow.
* [Basilisk](http://basilisk.fr/), C,  free software program for the solution of partial differential equations on adaptive Cartesian meshes.
* [OpenFOAM](https://www.openfoam.com/), C++, Open Source Field Operation and Manipulation (OpenFOAM) C++ libraries.
* [OpenFVM](https://openfvm.sourceforge.net/), C++, general CFD solver for the flow in complex 3D geometries.
* [OpenFlower](https://openflower.sourceforge.net/index2.html), C++, Large Eddy Simulations of turbulent flows.
* [FFR](https://www.eng.hokudai.ac.jp/labo/fluid/ssl/download/download.htm), Fortran
* [code_saturne](https://www.code-saturne.org/cms/web/), Fortran, developed primarily by EDF for CFD applications.
* [T-Flows](https://github.com/DelNov/T-Flows), Fortran, CFD program for simulation of turbulent, single and multiphase flows.
* [UCNS3D](https://github.com/ucns3d-team/UCNS3D), Fortran, versatile 2D and 3D unstructured CFD framework for a wide-range of compressible flow problems.
* [Clawpack](https://www.clawpack.org/), Fortran/Python, collection of finite volume methods for linear and nonlinear hyperbolic systems of conservation laws.

### FEM
* [nektar++](https://www.nektar.info/), C++, finite element package designed to allow one to construct efficient classical low polynomial order h-type solvers.
* [deal.II](https://dealii.org/), C++, solution of partial differential equations using adaptive finite elements.
* [Fluidity](https://fluidityproject.github.io/), C++, computational fluid dynamics code with adaptive unstructured mesh capabilities.
* [HiFlow](https://emcl-gitlab.iwr.uni-heidelberg.de/hiflow3.org/hiflow3/), C++, multi-purpose finite element software.
* [FrontFlow_blue](http://www.ciss.iis.u-tokyo.ac.jp/dl/download/permission.php?group_id=7), Fortran
* [FEATFLOW](https://wwwold.mathematik.tu-dortmund.de/~featflow/en/index.html), Fortran, incompressible flow in 2D and 3D.
* [elmer](https://www.elmerfem.org/blog/), Fortan/C, finite element software for multiphysical problems.
* [JuliaFEM](https://github.com/JuliaFEM/JuliaFEM.jl), Julia, open-source software for reliable, scalable, distributed Finite Element Method.
* [Firedrake](https://www.firedrakeproject.org/), Python·Cython, automated system for the solution of partial differential equations using the finite element method.
* [GetFEM](https://getfem.readthedocs.io/en/latest/), Python, Matlab, solving potentially coupled systems of linear and nonlinear partial differential equations.

### Compressible flow
* [TYCHO](https://www.tycho-cfd.at/), C, multidimensional (1D/2D and 3D) compressible hydrodynamics code.
* [HiFiLES](https://github.com/HiFiLES/HiFiLES-solver), C/C++, High Fidelity Large Eddy Simulation software package.
* [SU2](https://su2code.github.io/), C++, multiphysics analysis and design optimization software.
* 

function [n,bornes,ref,ldr,nr,ldi,ni]=retindice(ld,parm,method,varargin)
%  function [n,bornes,ref,ldr,nr,ldi,ni]=retindice(ld,parm,method);
%  indice de materiaux fonction de la longueur d'onde ld EN MICRONS
%  bornes:domaine de validite (quand interploation retourne nan à l'exterieur, sauf si les bornes sont entre parenthese s)
%  ref:references(chaine de characteres)
%  ldr,nr,ldi,ni:valeurs tabulees dans le cas de tabulations
%--------------------------------------------------------------------------------------------------------------------------
%   parm  | milieu  | bornes en lambda | references
%--------------------------------------------------------------------------------------------------------------------------
%    1    | argent  |  .2      10      | P.Lalanne:  Electronic Handbook of Optical Constant of solid version 1.0 HOC
%
%    1.28 |         | .000124  .04592  | JPH:   Electronic Handbook of Optical Constant of solid version 1.0 HOC
%    1.35 |         | .04592   .3647   | JPH:   Electronic Handbook of Optical Constant of solid version 1.0 HOC
%    1.33 |         | .3757    2.066   | JPH:   Electronic Handbook of Optical Constant of solid version 1.0 HOC
%    1.14 |         | 1.265    9.919   | JPH:   Electronic Handbook of Optical Constant of solid version 1.0 HOC
%    1.4  |         |  .75     .9      | TRES MAUVAIS ***  modele de Drude 'comment' (assez mediocre ...)
%    1.5  |         |  .75     .9      | modele de Drude benchmark
%    1.7  |         |  0.1879,1.9372   | P.B. Johnson. R.W. Christy. Phys. Rev. B6. 4370 (1972)
%    1.71 |         |  (        )      | modele de Drude donné par P.B. Johnson. R.W. Christy transmis par Christophe Sauvan
%    1.72 |         |  0.1879,2.0664   | Modele de Lorentz à 4 termes Normander 'Effective dielectric function for FDTD...' Chemical Physics Letters 446 (2007) 115-118 Transmis par JJG
%    1.8  |         |  (        )      | modele de Drude donné par A Archambault,F Marquier, JJ Greffet PRB 82 ,035411 (2010)
%    1.91 |         |  (.25,12.4)      | Modele de Lorentz-Drude A.D.Rakic,A B.Djurisic,J.M.Elazar,M.L.Majewski Appl.Opt.V37,5271-5283(1998)
%    1.92 |         |  (.25,12.4)      | Modele de Brendel_Bormann A.D.Rakic,A B.Djurisic,J.M.Elazar,M.L.Majewski Appl.Opt.V37,5271-5283(1998) 
%    1.6  |         |                  | modele de Drude Ordal transmis par P Ben-Abdallah
%    1.61 |         |                  |  correction JPH 2 2014 d'apres le  Palik pour grand ld
%--------------------------------------------------------------------------------------------------------------------------
%    2    |  or     | .15      15      | P.Lalanne:  Palik et Electronic Handbook of Optical Constant of solid version 1.0 HOC
%
%    2.01 |         | .5166    .9537   | JT:    M.C. Theye, Phys. Rev. B2, 3060 (1970)
%    2.02 |         | .521     .985    | obsolete voir 2.7     JT:    P.B. Johnson. R.W. Christy. Phys. Rev. B6. 4370 (1972)
%    2.03 |         | .5       .95     | JT:    K. Weiss. Z. Naturforscher 3a. 143 (1948)	
%    2.04 |         | .55      .95     | JT:    Corn   http://www.sopra   sa.com/index2.htm 
%    2.05 |         | .51662   .99191  | JT:    http://www.sopra   sa.com/index2.htm
%    2.06 |         | .4       .85     | Guillaume sat
%
%    2.28 |         | .000125  .0124   | JPH:   Electronic Handbook of Optical Constant of solid version 1.0 HOC
%    2.39 |         | .01409   .04768  | JPH:   Electronic Handbook of Optical Constant of solid version 1.0 HOC
%    2.10 |         | .047     .2      | JPH:   Electronic Handbook of Optical Constant of solid version 1.0 HOC
%    2.27 |         | .2066    2.066   | JPH:   Electronic Handbook of Optical Constant of solid version 1.0 HOC
%    2.14 |         | 1.265    4       | JPH:   Electronic Handbook of Optical Constant of solid version 1.0 HOC
%
%    2.2  |         | .20664   2.479   | H.Sauer:    ref='???????? ';
%    2.3  |         | .5       1       |  approximation polynomiale de degre 4 de Corn   (2.04)
% 
%    2.4  |         |   ?      ?       | TRES MAUVAIS ***  modele de Drude (mauvaise approximation)
%    2.5  |         |  (.5 1.1    ^)   | modele de Drude Lorentz: F.Kaminski,V.Sandoghdar,and M.Agio Journal of Computional and Theoretical Nanoscience Vol 4 635-643 2007
%
%    2.7  |         |  0.1879,1.9372   | P.B. Johnson. R.W. Christy. Phys. Rev. B6. 4370 (1972)
%    2.72 |         |  0.1879,2.0664   | Modele de Lorentz à 4 termes Normander 'Effective dielectric function for FDTD...' Chemical Physics Letters 446 (2007) 115-118 Transmis par JJG
%    2.91 |         |  (.25, 6.2)      | Modele de Lorentz-Drude A.D.Rakic,A B.Djurisic,J.M.Elazar,M.L.Majewski Appl.Opt.V37,5271-5283(1998)
%    2.92 |         |  (.25, 6.2)      | Modele de Brendel_Bormann A.D.Rakic,A B.Djurisic,J.M.Elazar,M.L.Majewski Appl.Opt.V37,5271-5283(1998) 
%    2.6  |         |                  |  TRES MAUVAIS ***  modele de Drude Ordal transmis par P Ben-Abdallah
%    2.61 |         |                  |  correction JPH 2 2014 d'apres le  Palik pour grand ld
%--------------------------------------------------------------------------------------------------------------------------
%    3    | Chrome  | 0.20664   1.239  | 
%   3.07  | Chrome  |  0.1879,1.9372   | P.B. Johnson. R.W. Christy. Phys. Rev. B6. 5059
%   20.7  |      idem
%    3.26 | Chrome  | 0.01   20 |Palick  et Electronic Handbook of Optical Constant of solid version 1.0 HOC transmis par Mathieu
%    3.7  |         | 0.0100     30.995| ????  Concatenation of: Cr_llnl_cxro + Cr_palik   http://www.luxpop.com????? pas compatible avec  3 et 20.7
%    3.91 |         |  (.25,  62)      | Modele de Lorentz-Drude A.D.Rakic,A B.Djurisic,J.M.Elazar,M.L.Majewski Appl.Opt.V37,5271-5283(1998)
%    3.92 |         |  (.25,  62)      | Modele de Brendel_Bormann A.D.Rakic,A B.Djurisic,J.M.Elazar,M.L.Majewski Appl.Opt.V37,5271-5283(1998) 
%--------------------------------------------------------------------------------------------------------------------------
%    4    | verre   |                  |  sans dispersion 
%    4.1  | SF11    |                  | 
%    4.2  | SF10    |                  | 
%    4.3  | BK7     |                  |
%    4.31 | SiO2    |                  | Sellmeier equation Wikipedia ajout 1 2013
%    4.4  | SiO2    | 3.7070, .2130    | extrait de la bibliothèque CodeV
%    4.41 | SiO2    | 1.5385, 20       | Palik transmis par Simon
%    4.42 | SiO2    |                  | modele de Lorentz d'apres 4.41 fité par JPH valable dans le domaine  1.5385, 20
%    4.45 |MoSiO2   |0.0001239   .41326|  Sébastien de Rossi. Calculé à partir des 'atomic scattering factors' de (6.31g/cm3) CXRO <http://www-cxro.lbl.gov> et LLNL <http://www-phys.llnl.gov/V_Div/scattering/asf.html>.
%    4.5  |Silicium |   1.1    5.6     | Frey Leviton Madison,Proc,SPIE Vol,6273,6272J'  
%    4.53 |Silicium |   .25    1.45    | Green,M.A.,and M.J.Keevers,Optical properties of intrinsic silicon at 300 K,Progress in Photovoltaics:Research and Applications,vol.3,issue 3,pp.189-192,1995 (transmis par Anthony Jouanin)
%    4.52 |Silicium |0.0001239   3332. |  Sébastien de Rossi      Si_llnl_cxro + Si_palik   
%    4.51 |Germanium|   1.9    5.5     | Frey Leviton Madison,Proc,SPIE Vol,6273,6272J'  
%--............ajout Stephane Collin
%    4.601|Al20Ga80As   |  0.5    2.48     | transmis par S.Collin from PGA {Properties of Gallium Arsenide} {S. Adachi} {Emis Datareviews Series No 7, London : INSPEC, 1993} - added by  NV (10/2012)
%   4.6021|Al30Ga70As   |  1.24e-4 100     | Palik transmis par S.Collin - added by  NV (11/2012)
%    4.602|Al315Ga685As |  0.5    2.48     | transmis par S.Collin from PGA {Properties of Gallium Arsenide} {S. Adachi} {Emis Datareviews Series No 7, London : INSPEC, 1993} - added by  NV (10/2012)
%    4.603|Al50Ga50As   |  0.5    2.48     | transmis par S.Collin from PGA {Properties of Gallium Arsenide} {S. Adachi} {Emis Datareviews Series No 7, London : INSPEC, 1993} - added by  NV (10/2012)
%    4.604|Al70Ga70As   |  0.5    2.48     | transmis par S.Collin from PGA {Properties of Gallium Arsenide} {S. Adachi} {Emis Datareviews Series No 7, London : INSPEC, 1993} - added by  NV (10/2012)
%    4.605|Al80Ga20As   |  0.21   2.48     | recollï¿½ (NV) de deux tables de Palik - added by NV
%    4.606|AlAs         |  0.5    2.48     | transmis par S.Collin from PGA {Properties of Gallium Arsenide} {S. Adachi} {Emis Datareviews Series No 7, London : INSPEC, 1993} - added by  NV (10/2012)
%    4.701|GaAs         |  0.008  1000     | Palik transmis par S.Collin -    added by  NV (10/2012)
%    4.91 |GaAs:Si      |  .3     3.0036   | n+ Silicon doped GaAs 2*10^18    J. Electromagnetic Analysis & Applications, 2010, 2, 357-361
%    4.13 |GaSb         |  0.2067 100      | cubic - Palik & SOPRA - re-added by NV (10/2012)
%--............fin ajout Stephane Collin
%    4.6  |AlxGA1-xAs|                 | Optoélectronique, Emmanuel Rosencher et Borge Vinter, Masson, chapitre 12.7, P332'; 
%    4.61 |PxAs1-xGa|                  | Optoélectronique, Emmanuel Rosencher et Borge Vinter, Masson, chapitre 12.7, P332'; 
%    4.62 |GaxIn1-xP|                  | Optoélectronique, Emmanuel Rosencher et Borge Vinter, Masson, chapitre 12.7, P332'; 
%    4.63 |GaAs dopé|                  | formule de  Simon Vassant Palik  Electronic Handbook of Optical Constant of solid ISBN 0-12-544420-60 p429
%    4.71 |Saphir no| 1.1098 1.6932    | Laurent FREY Handbook of Optical   Constants of Solids III, Ed. D. Palik (1998) transmis par P Lalanne
%    4.72 |Saphir ne|                  |
%    4.7  |GaAs     |.08856      17.46 | indice de  GaAs Palik p 434 440 4.701 est preferable car couvre un domaine plus grand en completant
%    4.8  |TiO2     |.35            .8 | 'Yamada Deposition at low substrate temperatures of hight quality TiO2 Films by radical beam assistef evaporation Appl Opt 1999 38 6638-6641 courbe Fig2 300K fitee par JPH ';
%   4.81  |TiO2     | 0.001      4.5   |Palick  et Electronic Handbook of Optical Constant of solid version 1.0 HOC, transmis par mathieu
%   4.82  |TiO2     |                  |Sellmeier equation transmis par  Mathieu http://refractiveindex.info
%   4.9   |MgF2     |                  |Sellmeier equation transmis par Mathieu  http://refractiveindex.info   
%  4.11   |Al2O3    |                  |Sellmeier equation transmis par Mathieu  http://refractiveindex.info   
%  4.14   |SiO2     |                  |Sellmeier equation transmis par Mathieu  http://refractiveindex.info
%  4.12   |ZnO      |                  |Sellmeier equation transmis par Mathieu  http://refractiveindex.info
%  4.15   | Si3N4   |0.18    4.8       | Palik (?) transmis par Mathieu';
%--------------------------------------------------------------------------------------------------------------------------
%    5    | eau     |                  | 
%    5.1  | eau SD  |                  |  sans dispersion 
%--------------------------------------------------------------------------------------------------------------------------
%    6    | air     |                  |  sans dispersion 
%    7    | bio     |                  |  sans dispersion 
%    8    | résine  |                  |  sans dispersion 
%    9    |polypyrrole|                |  sans dispersion 
%--------------------------------------------------------------------------------------------------------------------------
%    10   |aluminium| 0.00415  30.99   | http://cvs.savannah.nongnu.org/viewvc/freesnell/al.nk?root=freesnell&view=markup  project freesnell - View of /freesnell/al.nk 23 Octobre 2004 (Transmis par Inbal Friedler et P.Lalanne)
%    10.1 |         |                  | TRES MAUVAIS *** modele de Drude PhysRevB.68.245405
%   10.91 |         |  (        )      | Modele de Lorentz-Drude A.D.Rakic,A B.Djurisic,J.M.Elazar,M.L.Majewski Appl.Opt.V37,5271-5283(1998)
%   10.92 |         |  (        )      | Modele de Brendel_Bormann A.D.Rakic,A B.Djurisic,J.M.Elazar,M.L.Majewski Appl.Opt.V37,5271-5283(1998) 
%   10.6  |         |                  |  TRES MAUVAIS *** modele de Drude Ordal transmis par P Ben-Abdallah 
%--------------------------------------------------------------------------------------------------------------------------
%    12   |  cuivre |.000137756,9.53692| |  A VERIFIER | Handbook of Optical Constants of Solids, Ed. by Edward D. Palik,http://www.luxpop.com
%    12.7 |         |  0.1879,1.9372   | P.B. Johnson. R.W. Christy. Phys. Rev. B6. 4370 (1972)
%   12.91 |         |  (        )      | Modele de Lorentz-Drude A.D.Rakic,A B.Djurisic,J.M.Elazar,M.L.Majewski Appl.Opt.V37,5271-5283(1998)
%   12.92 |         |  (        )      | Modele de Brendel_Bormann A.D.Rakic,A B.Djurisic,J.M.Elazar,M.L.Majewski Appl.Opt.V37,5271-5283(1998) 
%   12.6  |         |  ( 55   ..)      |   TRES MAUVAIS *** modele de Drude Ordal transmis par P Ben-Abdallah 
%   12.61 |         |                  |  correction JPH 2 2014
%--------------------------------------------------------------------------------------------------------------------------
%    13   | Tungsten|0.000012398,24.7960|  'W_llnl_cxro + W_palik'; http://www.luxpop.com (corrigé en 11 2013 )
%  13.24  |Tungsten | 0.1       24.80  | Electronic Handbook of Optical Constant of solid version 1.0 HOC transmis par Mathieu
%  13.241 |         | 0.1       4.5    |Electronic Handbook of Optical Constant of solid version 1.0 HOC transmis par Mathieu
%   13.91 |         |  (.25, 12.4)     | Modele de Lorentz-Drude A.D.Rakic,A B.Djurisic,J.M.Elazar,M.L.Majewski Appl.Opt.V37,5271-5283(1998)
%   13.92 |         |  (.25, 12.4)     | Modele de Brendel_Bormann A.D.Rakic,A B.Djurisic,J.M.Elazar,M.L.Majewski Appl.Opt.V37,5271-5283(1998) 
%--------------------------------------------------------------------------------------------------------------------------
%    14.7 | Fer     |  0.1879,1.9372   | P.B. Johnson. R.W. Christy. Phys. Rev. B6. 5059
%--------------------------------------------------------------------------------------------------------------------------
%    15.7 | Cobalt  |  0.1879,1.9372   | P.B. Johnson. R.W. Christy. Phys. Rev. B6. 5059
%    15.6 |         |                  |  TRES MAUVAIS ***  modele de Drude Ordal transmis par P Ben-Abdallah
%--------------------------------------------------------------------------------------------------------------------------
%    16.7 | Nickel  |  0.1879,1.9372   | P.B. Johnson. R.W. Christy. Phys. Rev. B6. 5059
%   16.91 |         |  (.25, 6.2)      | Modele de Lorentz-Drude A.D.Rakic,A B.Djurisic,J.M.Elazar,M.L.Majewski Appl.Opt.V37,5271-5283(1998)
%   16.92 |         |  (.25, 6.2)      | Modele de Brendel_Bormann A.D.Rakic,A B.Djurisic,J.M.Elazar,M.L.Majewski Appl.Opt.V37,5271-5283(1998) 
%--------------------------------------------------------------------------------------------------------------------------
%    17.7 |Palladium|  0.1879,1.9372   | P.B. Johnson. R.W. Christy. Phys. Rev. B6. 5059
%   17.91 |         |  (.25,12.4)      | Modele de Lorentz-Drude A.D.Rakic,A B.Djurisic,J.M.Elazar,M.L.Majewski Appl.Opt.V37,5271-5283(1998)
%   17.92 |         |  (.25,12.4)      | Modele de Brendel_Bormann A.D.Rakic,A B.Djurisic,J.M.Elazar,M.L.Majewski Appl.Opt.V37,5271-5283(1998) 
%--------------------------------------------------------------------------------------------------------------------------
%    18   |  Titane | 0.000123  200    | Electronic HandBook of Optical Constants of Solids (HOC) transmis par JC Rodier
%    11   |      idem  
%    18.7 |         |  0.1879,1.9372   | P.B. Johnson. R.W. Christy. Phys. Rev. B6. 5059
%   18.91 |         |  (.25,  30)      | Modele de Lorentz-Drude A.D.Rakic,A B.Djurisic,J.M.Elazar,M.L.Majewski Appl.Opt.V37,5271-5283(1998)
%   18.92 |         |  (.25,  30)      | Modele de Brendel_Bormann A.D.Rakic,A B.Djurisic,J.M.Elazar,M.L.Majewski Appl.Opt.V37,5271-5283(1998) 
%   18.6  |         |  (13,   ...)     |  TRES MAUVAIS ***  modele de Drude Ordal transmis par P Ben-Abdallah 
%   11.6  |      idem   
%   18.61 |         |                  |  correction JPH 2 2014
%--------------------------------------------------------------------------------------------------------------------------
%    19.7 | Vanadium|  0.1879,1.9372   | P.B. Johnson. R.W. Christy. Phys. Rev. B6. 5059
%--------------------------------------------------------------------------------------------------------------------------
%    21.7 |Manganese|  0.1879,1.9372   | P.B. Johnson. R.W. Christy. Phys. Rev. B6. 5059
%--------------------------------------------------------------------------------------------------------------------------
%    22   |MgO      |  12.2, 43.4      | Handbook of Optical Constants of Solids, Ed. by Edward D. Palik, transmis par Mondher
%    22.1 |MgO      | 0.01653   20     | Palick  et Electronic Handbook of Optical Constant of solid
%--------------------------------------------------------------------------------------------------------------------------
%    23   |SiC      |                  | modele de Lorentz donnees Palik transmis par:J le Gall M Olivier JJ Greffet Phys Rev B vol 55 number 15 15 april 1997
%    23.2 | SiC     |                  | Modèle de Lorentz pour SiC  non poli transmis par Mathieu
%    23.3 | SiC     |                  | Modèle de Lorentz pour SiC poli    transmis par Mathieu
%    23.1 |SiC      |0.04133,   25     | Palik hexagonal SiC, ordinary ray.  L'accord avec 23 est remarquable dans le domaine 5 25 mais Palik manque de points dans la zone 12.4 12.8: 23 est preferable 
%    23.12| SiC     | 0.262  , 10      | Mesures  ellipso transmis par Mathieu
% 	 23.11| SiC     | 0.262  , 10      | Fit Mesures      transmis par Mathieu
%--------------------------------------------------------------------------------------------------------------------------
%   24    |Platine  | 0.3757  12.4     | Palik transmis par Daniele Costantini
%   24.91 |         |  (.25, 12.4)     | Modele de Lorentz-Drude A.D.Rakic,A B.Djurisic,J.M.Elazar,M.L.Majewski Appl.Opt.V37,5271-5283(1998)
%   24.92 |         |  (.25, 12.4)     | Modele de Brendel_Bormann A.D.Rakic,A B.Djurisic,J.M.Elazar,M.L.Majewski Appl.Opt.V37,5271-5283(1998) 
%--------------------------------------------------------------------------------------------------------------------------
%   25.91 |Beryllium|  (.25,  62)      | Modele de Lorentz-Drude A.D.Rakic,A B.Djurisic,J.M.Elazar,M.L.Majewski Appl.Opt.V37,5271-5283(1998)
%   25.92 |         |  (.25,  62)      | Modele de Brendel_Bormann A.D.Rakic,A B.Djurisic,J.M.Elazar,M.L.Majewski Appl.Opt.V37,5271-5283(1998)
%--------------------------------------------------------------------------------------------------------------------------
%    27   | Tantale | 0.01      125    | Palik, transmis par Mathieu
%--------------------------------------------------------------------------------------------------------------------------
%    26   |Molybdene| 0.01      12.40  |Electronic Handbook of Optical Constant of solid version 1.0 HOC, transmis par Mathieu
%    26.1 |         | 0.01      4.5    |Electronic Handbook of Optical Constant of solid version 1.0 HOC, transmis par Mathieu
%
%-------------------------------------------------------------------------
% --...........ajout Stephane Collin    modif par rapport à Stephane 23 -> 28 car 23 est SiC
%                            AR Coatings 
%    28.1   | c-ZnS                    | 0.0001240 666.7   | Palik - added by NV (03/2012)
%    28.2   | Si3N4                    | 0.05166   1.24    | polycristalline - Palik - added by NV (03/2012)
%    28.31  | MgF2 ordinary index      |   0.0459  2000    | tetragonal - Palik - added by NV (03/2012)
%    28.32  | MgF2 extraordinary index |   0.0459  2000    | tetragonal - Palik - added by NV (03/2012)
%    28.41  | TiO2 ordinary index      |   0.103   13.2    | tetragonal (rutile) - Palik - added by NV (10/2012)
%    28.42  | TiO2 extraordinary index |   0.103   13.2    | tetragonal (rutile) - Palik - added by NV (10/2012)
%    28.43  | TiO2 SOPRA 1             |   0.18    1.5     | tetragonal (rutile) - SOPRA - added by NV (10/2012)
%    28.44  | TiO2 SOPRA 2             |   0.22    2.42    | tetragonal (rutile) - SOPRA - added by NV (10/2012)
%-- INTERPOLES 
%   101.7  |             |  .3     1.91     | S.Collin 2_2010
%   102.7  |             |  .3     1.91     | S.Collin  2_2010
%   103.7  |             |  .3     1.91     | S.Collin  2_2010
%   104.7  |GaAs         |  .3     1.91     | S.Collin  2_2010
%   104.8  |GaAs:C       |  0.3   1.795     | p+ Carbon doped GaAs  2*10^19
%--....................
%    114   | ZnO         | 0.3    1.045     | Palik 
%    114.2 | ZnO i       | 0.355  , 1.9972  | Ellipsométrie Xavier Lafosse LPN
%    114.1 | Zno:Al      | 0.1     1.91     | S.Collin transmis par Ch.Sauvan 2_2010
%    114.3 | ZnO:Al      | 0.355 , 1.9972   | Ellipsométrie Xavier Lafosse LPN
%--....................
%    115   | CdS (on glass)          |  0.3   1.91   | Palik
%    115.1 | CdS (thin film on CIGS) |  0.35  1.39   | Malmstrom, cas particul
%--....................
%    116   | Mo          |  0.3    3.0036 | Palik
%--....................
%   117     | CIGS23    |  0.355  1.455   | Orgassa via Jonas Malmstrom
%   117.1   | CIGS31    |  0.353  1.455   | Paulson via Jonas Malmstrom
%--....................
%   118     | microSi   |0.353, 0.83      | LPICM
%   118.1   | a-Si      |0.355, 0.85      | LPICM
%   118.2   | a-Si i    |0.3       1.47   | LPICM (01/2012)
%   118.3   | a-SiC p   |0.3       1.47   | LPICM (01/2012)
%--....................
%   119     | Custom Alu/Argent  |  0.3  2    | Xavier Lafosse LPN Ellipso
%--....................
%   120     | InP       |  .3  1  | Palik, Cubic Indium Phosphide
%--....................
%   121     | ITO       | 0.3  1.13 | LPICM (01/2012)
%   122     | ITO       | 0.3  1.12 | LPICM (02/2012)
% --.......fin ajout Stephane Collin
%-------------------------------------------------------------------------
% --...........ajout Mathieu
%--------------------------------------------------------------------------------------------------------------------------
%   30  |   Glass_HEF | 0.3     2   | Mesure R & T, fit avec SCOUT
%   31  |   W HEF   |   0.3     2   | Mesures R & T, fit avec SCOUT
%   32  |   MgO HEF |   0.3     2   | Mesures R & T, fit avec SCOUT
%   33  |   TiO2 HEF|   0.3     2   | Mesures R & T, fit avec SCOUT
%--------------------------------------------------------------------------------------------------------------------------
%   40  |   MgO_MD  | 0.1     22   | MgO Simulation material design
%   41  |   SiC_MD  | 0.2     22   | SiC cubic Simulation material design
%   42  |   TiO2_MD | 0.2     22   | TiO2 rutile Simulation material design
%   43  |   W_MD    | 0.2     6    | W Simulation material design
%   44  |   WO2_MD  | 0.2     11    | WO2 Simulation material design
% --...........fin ajout Mathieu
%
%
% Par Liste Alphabetique: voir retindice_liste_comparative dans rettest nv_retindice
% -----------------------------------------------------------------------------------------------------
% | Materiel            |   parm        |BORNES( Inf =analytique)   | references
% -----------------------------------------------------------------------------------------------------
% | AL                  | 10            | 0.00415      30.99        |  Savannah. 
% | Ag                  | 1             | 0.2          10           |   Electronic Handbook of Optical Const 
% | Ag                  | 1.28          | 0.000124     0.04592      |    JPH:   Electronic Handbook of Optic 
% | Ag                  | 1.35          | 0.04592      0.3647       |    JPH:   Electronic Handbook of Optic 
% | Ag                  | 1.33          | 0.3757       2.066        |    JPH:   Electronic Handbook of Optic 
% | Ag                  | 1.14          | 1.265        9.919        |    JPH:   Electronic Handbook of Optic 
% | Ag                  | 1.4           | 0.75         0.9          |   ASSEZ MEDIOCRE *** modele de Drude ' 
% | Ag                  | 1.5           | 0.4          2            |   TRES MAUVAIS *** modele de Drude 'be 
% | Ag                  | 1.7           | 0.18785      1.9372       |   P.B. Johnson. R.W. Christy. Phys. Re 
% | Ag                  | 1.71          | -Inf         Inf          |   modéle de Drude donné par P.B. Johns 
% | Ag                  | 1.72          | -Inf         Inf          |  Modele de Lorentz à 4 termes Normande 
% | Ag                  | 1.8           | -Inf         Inf          |  modele de Drude donné par A Archambau 
% | Ag                  | 1.91          | -Inf         Inf          |  A.D.Rakic,A B.Djurisic,J.M.Elazar,M.L 
% | Ag                  | 1.92          | -Inf         Inf          |  A.D.Rakic,A B.Djurisic,J.M.Elazar,M.L 
% | Ag                  | 1.6           | -Inf         Inf          |  TRES MAUVAIS ***modele de Drude Ordal 
% | Ag                  | 1.61          | -Inf         Inf          |  modele de Drude Ordal transmis par P  
% | Ag                  | 101.7         | 0.3          1.91         |   Donnees interpollees par S Collin co 
% | Air                 | 6             | -Inf         Inf          |  sans dispersion 
% | Al                  | 10.1          | -Inf         Inf          |  modele de Drude PhysRevB.68.245405 
% | Al                  | 10.91         | -Inf         Inf          |  A.D.Rakic,A B.Djurisic,J.M.Elazar,M.L 
% | Al                  | 10.92         | -Inf         Inf          |  A.D.Rakic,A B.Djurisic,J.M.Elazar,M.L 
% | Al                  | 10.6          | -Inf         Inf          |  TRES MAUVAIS ***modele de Drude Ordal 
% | Al20Ga80As          | 4.601         | 0.49633      10           |  : S. Adachi : Properties of G 
% | Al2O3               | 4.11          | -Inf         Inf          |  Sellmeier equation, transmis par M 
% | Al2O3               | 4.12          | -Inf         Inf          |  Sellmeier equation, transmis par M 
% | Al30Ga70As          | 4.6021        | 0.000124     100          |   Palik optical constants : Al 
% | Al315Ga685As        | 4.602         | 0.49633      2.4816       |  : S. Adachi : Properties of 
% | Al50Ga50As          | 4.603         | 0.49633      2.4816       |  : S. Adachi : Properties of G 
% | Al70Ga30As          | 4.604         | 0.49633      2.4816       |  : S. Adachi : Properties of G 
% | Al80Ga20As          | 4.605         | 0.2066       2.48         |  Palik Optical Constants (deux 
% | AlAs                | 4.606         | 0.49633      2.4816       |  : S. Adachi : Properties of Gallium 
% | AlxGA1-xAs          | 4.6           | -Inf         Inf          |  Optoélectronique, Emmanuel Ro 
% | Au                  | 2             | 0.15         15           |    Palik  et Electronic Handbook of Op 
% | Au                  | 2.01          | 0.5166       0.9537       |    JT:  M.C. Theye, Phys. Rev. B2, 306 
% | Au                  | 2.02          | 0.521        0.985        |    JT:     P.B. Johnson. R.W. Christy. 
% | Au                  | 2.03          | 0.5          0.95         |   JT:     K. Weiss. Z. Naturforscher 3 
% | Au                  | 2.04          | 0.55         0.94         |    JT:     Corn	 http://corndog.chem.w 
% | Au                  | 2.05          | 0.51662      0.99191      |   JT:   http://www.sopra   sa.com/inde 
% | Au                  | 2.06          | 0.4          0.85         |   guillaume sat 
% | Au                  | 2.28          | 0.000125     0.0124       |    JPH:   Electronic Handbook of Optic 
% | Au                  | 2.39          | 0.01409      0.04768      |    JPH:   Electronic Handbook of Optic 
% | Au                  | 2.1           | 0.047        0.2          |    JPH:   Electronic Handbook of Optic 
% | Au                  | 2.27          | 0.2066       2.066        |    JPH:   Electronic Handbook of Optic 
% | Au                  | 2.14          | 1.265        4            |    JPH:   Electronic Handbook of Optic 
% | Au                  | 2.2           | 0.20664      2.479        |  ????????  
% | Au                  | 2.3           | 0.55         0.95         |   approximation polynomiale de degre 4 
% | Au                  | 2.4           | -Inf         Inf          |  TRES MAUVAIS *** modele de Drude  ben 
% | Au                  | 2.5           | -Inf         Inf          |  modele de Drude Lorentz: F.Kaminski,V 
% | Au                  | 2.7           | 0.18785      1.9372       |   P.B. Johnson. R.W. Christy. Phys. Re 
% | Au                  | 2.72          | -Inf         Inf          |  Modele de Lorentz à 4 termes Normande 
% | Au                  | 2.91          | -Inf         Inf          |  A.D.Rakic,A B.Djurisic,J.M.Elazar,M.L 
% | Au                  | 2.92          | -Inf         Inf          |  A.D.Rakic,A B.Djurisic,J.M.Elazar,M.L 
% | Au                  | 2.6           | -Inf         Inf          |  TRES MAUVAIS ***modele de Drude Ordal 
% | Au                  | 2.61          | -Inf         Inf          |  modele de Drude Ordal transmis par P  
% | Au                  | 102.7         | 0.3          1.91         |   Donnees interpollees par S Collin co 
% | BK7                 | 4.3           | -Inf         Inf          |  bornes a preciser Sellmeier equation 
% | Be                  | 25.91         | -Inf         Inf          |  A.D.Rakic,A B.Djurisic,J.M.Elazar,M.L 
% | Be                  | 25.92         | -Inf         Inf          |  A.D.Rakic,A B.Djurisic,J.M.Elazar,M.L 
% | Bio                 | 7             | -Inf         Inf          |  sans dispersion 
% | CIGS23              | 117           | 0.355        1.455        |   Donnees interpollees par S Colli 
% | CIGS31              | 117.1         | 0.355        1.53         |   Donnees interpollees par S Colli 
% | CR                  | 3.7           | 0.0099795    30.995       |  Concatenation of: Cr_llnl_cxro + Cr_p 
% | CdS                 | 115           | 0.3          1.91         |   Donnees interpollees par S Collin c 
% | CdS2                | 115.1         | 0.355        1.39         |   Donnees interpollees par S Collin  
% | Co                  | 15.7          | 0.18785      1.9372       |    P.B. Johnson. R.W. Christy. Phys. R 
% | Co                  | 15.6          | -Inf         Inf          |   TRES MAUVAIS ***modele de Drude Orda 
% | Cr                  | 3             | 0.20664      1.239        |   ????????  
% | Cr                  | 3.07          | 0.18785      1.9372       |   P.B. Johnson. R.W. Christy. Phys. Re 
% | Cr                  | 20.7          | 0.18785      1.9372       |   P.B. Johnson. R.W. Christy. Phys. Re 
% | Cr                  | 3.26          | 0.01         20           |    Palick  et Electronic Handbook of O 
% | Cr                  | 3.91          | -Inf         Inf          |  A.D.Rakic,A B.Djurisic,J.M.Elazar,M.L 
% | Cr                  | 3.92          | -Inf         Inf          |  A.D.Rakic,A B.Djurisic,J.M.Elazar,M.L 
% | Cr                  | 103.7         | 0.3          1.91         |   Donnees interpollees par S Collin co 
% | Cu                  | 12            | 0.00013776   9.5369       |  Handbook of Optical Constants of Soli 
% | Cu                  | 12.7          | 0.18785      1.9372       |   P.B. Johnson. R.W. Christy. Phys. Re 
% | Cu                  | 12.91         | -Inf         Inf          |  A.D.Rakic,A B.Djurisic,J.M.Elazar,M.L 
% | Cu                  | 12.92         | -Inf         Inf          |  A.D.Rakic,A B.Djurisic,J.M.Elazar,M.L 
% | Cu                  | 12.6          | -Inf         Inf          |  TRES MAUVAIS ***modele de Drude Ordal 
% | Cu                  | 12.61         | -Inf         Inf          |  modele de Drude Ordal transmis par P  
% | CustomAl            | 119           | 0.3          1.9972       |   Donnees interpollees par S Col 
% | Fe                  | 14.7          | 0.18785      1.9372       |   P.B. Johnson. R.W. Christy. Phys. Re 
% | GaAs                | 4.701         | 0.007999     1000         |  Palik Optical Constants 
% | GaAs                | 4.63          | -Inf         Inf          |   formule de  Simon Vassant Palik  E 
% | GaAs                | 4.7           | 0.08856      17.46        |    Palik p 434 440 
% | GaAs                | 104.7         | 0.3          1.91         |   Donnees interpollees par S Collin  
% | GaAs_C              | 104.8         | 0.3          1.795        |   Donnees interpollees par S Colli 
% | GaAs_dope_Si        | 4.91          | 0.3          3.0036       |  (n-type) J. Electromagnetic 
% | GaSb                | 4.13          | 0.2067       100          |  Palik Optical Constants 
% | GaxIn1-xP           | 4.62          | -Inf         Inf          |  Optoélectronique, Emmanuel Ros 
% | Ge                  | 4.5           | 1.1          5.6          |  Temperature dependent refractive inde 
% | Ge                  | 4.51          | 1.9          5.5          |  Temperature dependent refractive inde 
% | H2O                 | 5             | -Inf         Inf          |  avec dispersion, bornes a preciser,s 
% | H2O                 | 5.1           | -Inf         Inf          |  sans dispersion  
% | ITO                 | 121           | 0.3          1.13         |   Donnees interpollees par S Collin c 
% | ITO_2               | 122           | 0.3          1.12         |   Donnees interpollees par S Collin 
% | InP                 | 120           | 0.3          1            |   Donnees interpollees par S Collin c 
% | MgF2                | 4.9           | -Inf         Inf          |  Sellmeier equation, transmis par Ma 
% | MgF2_ne             | 28.32         | 0.0459       2000         |  Palik Optical Constants : MgF2 n 
% | MgF2_no             | 28.31         | 0.0459       2000         |  Palik Optical Constants : MgF2 n 
% | MgO                 | 22            | 12.2         43.47        |  Palik transmis par Mondher 
% | MgO                 | 22.1          | 0.01653      20           |    Palick  et Electronic Handbook of  
% | MgO                 | 32            | 0.3          2            |  HEF, mesures R&T, fit avec SCOUT 
% | MgO                 | 40            | 0.1          22           |  Material Design 
% | Mn                  | 21.7          | 0.18785      1.9372       |   P.B. Johnson. R.W. Christy. Phys. Re 
% | Mo                  | 26            | 0.001        12.4         |    Palick  et Electronic Handbook of O 
% | Mo                  | 26.1          | 0.001        4.5          |    Palick  et Electronic Handbook of O 
% | Mo                  | 116           | 0.3          3.0036       |   Donnees interpollees par S Collin co 
% | MoSiO2              | 4.45          | 0.00012398   0.41326      |  Sébastien de Rossi  CXRO <http:// 
% | Ni                  | 16.7          | 0.18785      1.9372       |   P.B. Johnson. R.W. Christy. Phys. Re 
% | Ni                  | 16.91         | -Inf         Inf          |  A.D.Rakic,A B.Djurisic,J.M.Elazar,M.L 
% | Ni                  | 16.92         | -Inf         Inf          |  A.D.Rakic,A B.Djurisic,J.M.Elazar,M.L 
% | Pd                  | 17.7          | 0.18785      1.9372       |   P.B. Johnson. R.W. Christy. Phys. Re 
% | Pd                  | 17.91         | -Inf         Inf          |  A.D.Rakic,A B.Djurisic,J.M.Elazar,M.L 
% | Pd                  | 17.92         | -Inf         Inf          |  A.D.Rakic,A B.Djurisic,J.M.Elazar,M.L 
% | Polypyrrole         | 9             | -Inf         Inf          |  sans dispersion 
% | Pt                  | 24            | 0.3757       43.47        |   Palik transmis par Daniele 
% | Pt                  | 24.91         | -Inf         Inf          |  A.D.Rakic,A B.Djurisic,J.M.Elazar,M.L 
% | Pt                  | 24.92         | -Inf         Inf          |  A.D.Rakic,A B.Djurisic,J.M.Elazar,M.L 
% | PxAs1-xGa           | 4.61          | -Inf         Inf          |  Optoélectronique, Emmanuel Ros 
% | Resine              | 8             | -Inf         Inf          |  sans dispersion 
% | SF10                | 4.2           | -Inf         Inf          |  bornes a preciser 
% | SF11                | 4.1           | -Inf         Inf          |  bornes a preciser 
% | Saphir_ne           | 4.72          | 1.014        1.6932       |  Handbook of Optical Constants  
% | Saphir_no           | 4.71          | 1.014        1.6932       |  Handbook of Optical Constants  
% | Si                  | 4.53          | 0.25         1.45         |   Green,M.A.,and M.J.Keevers,Optical p 
% | Si                  | 4.52          | 0.00012398   3332         |  llnl cxro + Si palik 
% | Si3N4               | 4.15          | 0.2          6            |  (Palik ? Mathieu) 
% | Si3N4               | 28.2          | 0.05166      1.24         |   Palik Optical constants I 
% | SiC                 | 23            | -Inf         Inf          |   modele de Lorentz donnees Palik tra 
% | SiC                 | 23.1          | 0.04133      25           |  Palik hexagonal SiC, ordinary ray ht 
% | SiC                 | 23.12         | 0.262        10           |  Mesure ellipso 
% | SiC                 | 23.11         | 0.262        10           |  Mesure ellipso fitée 
% | SiC                 | 41            | 0.2          13           |  cubic Material Design 
% | SiC_non_poli        | 23.3          | -Inf         Inf          |  Double oscillateur de Loren 
% | SiC_np              | 23.2          | -Inf         Inf          |  Double oscillateur de Lorentz 
% | SiO2                | 4.31          | -Inf         Inf          |   Silice_fondue bornes a preciser Se 
% | SiO2                | 4.4           | 3.707        0.213        |  avec dispersion, (extrait de la bib 
% | SiO2                | 4.41          | 1.5385       20           |  Palik transmis par Simon 
% | SiO2                | 4.42          | -Inf         Inf          |  modele de Lorentz d'apres 4.41 fité 
% | SiO2                | 4.14          | -Inf         Inf          |  Sellmeier equation,handbook of Opti 
% | Ta                  | 27            | 0.000124     125          |    Palick  et Electronic Handbook of O 
% | Ti                  | 18            | 0.000123     200          |  Electronic HandBook of Optical Consta 
% | Ti                  | 11            | 0.000123     200          |  Electronic HandBook of Optical Consta 
% | Ti                  | 18.7          | 0.18785      1.9372       |    P.B. Johnson. R.W. Christy. Phys. R 
% | Ti                  | 18.91         | -Inf         Inf          |  A.D.Rakic,A B.Djurisic,J.M.Elazar,M.L 
% | Ti                  | 18.92         | -Inf         Inf          |  A.D.Rakic,A B.Djurisic,J.M.Elazar,M.L 
% | Ti                  | 18.6          | -Inf         Inf          |   TRES MAUVAIS ***modele de Drude Orda 
% | Ti                  | 11.6          | -Inf         Inf          |   TRES MAUVAIS ***modele de Drude Orda 
% | Ti                  | 18.61         | -Inf         Inf          |   modele de Drude Ordal transmis par P 
% | TiO2                | 4.8           | 0.35         0.8          |  Yamada Deposition at low substrate  
% | TiO2                | 4.81          | 0.001        4.5          |    Palick  et Electronic Handbook of 
% | TiO2                | 4.82          | -Inf         Inf          |  Sellmeier equation, transmis par Ma 
% | TiO2                | 28.43         | 0.18         1.5          |   SOPRA  1 
% | TiO2                | 28.44         | 0.22         2.42         |   SOPRA  2 
% | TiO2                | 33            | 0.3          2            |  HEF, mesures R&T, fit avec SCOUT 
% | TiO2                | 42            | 0.2          13           |  rutile Material Design 
% | TiO2_extraordinary  | 28.42         | 0.103        13.2         |  Palik optical constan 
% | TiO2_ordinary       | 28.41         | 0.103        13.2         |  Palik optical constants :  
% | V                   | 19.7          | 0.18785      1.9372       |    P.B. Johnson. R.W. Christy. Phys. Re 
% | Verre               | 4             | -Inf         Inf          |  sans dispersion  
% | Verre               | 30            | 0.3          2            |  HEF, mesures R&T, fit avec SCOUT 
% | W                   | 13            | 1.2398e-005  24.796       |  W_llnl_cxro + W_palik 
% | W                   | 13.24         | 0.01         24.8         |    Palick  et Electronic Handbook of Op 
% | W                   | 13.241        | 0.01         4.5          |    Palick  et Electronic Handbook of Op 
% | W                   | 13.91         | -Inf         Inf          |  A.D.Rakic,A B.Djurisic,J.M.Elazar,M.L. 
% | W                   | 13.92         | -Inf         Inf          |  A.D.Rakic,A B.Djurisic,J.M.Elazar,M.L. 
% | W                   | 31            | 0.3          2            |  HEF, mesures R&T, fit avec SCOUT 
% | W                   | 43            | 0.2          6            |  Material Design 
% | WO2                 | 44            | 0.2          6            |  Material Design 
% | ZnO                 | 114           | 0.3          1.045        |   Donnees interpollees par S Collin c 
% | ZnO:Al2             | 114.3         | 0.355        1.9972       |   Donnees interpollees par S Coll 
% | ZnO_Al              | 114.1         | 0.3          1.91         |   Donnees interpollees par S Colli 
% | ZnOi2               | 114.2         | 0.355        1.9972       |   Donnees interpollees par S Collin 
% | ZnS_cubic           | 28.1          | 0.000124     666.7        |  Palik Optical constants I 
% | a_Si                | 118.1         | 0.355        0.825        |   Donnees interpollees par S Collin  
% | a_SiC_p             | 118.3         | 0.3          1.47         |   Donnees interpollees par S Coll 
% | a_Si_i              | 118.2         | 0.3          1.47         |   Donnees interpollees par S Colli 
% | microSi             | 118           | 0.355        0.83         |   Donnees interpollees par S Coll 
% ---------------------------------------------------------------------------------------------------
%
%                         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         %%%   LISTE SIMPLIFIEE    %%%
%                         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  function [n,bornes,ref,ldr,nr,ldi,ni]=retindice(ld,parm,method);
%  indice de materiaux fonction de la longueur d'onde ld EN MICRONS
%  bornes:domaine de validite (quand interploation retourne nan à l'exterieur, sauf si les bornes sont entre parenthese s)
%  ref:references(chaine de characteres)
%  ldr,nr,ldi,ni:valeurs tabulees dans le cas de tabulations
%--------------------------------------------------------------------------------------------------------------------------
%   parm  | milieu  | bornes en lambda | references
%--------------------------------------------------------------------------------------------------------------------------
%    1    | argent  |  .2      10      | P.Lalanne:  Electronic Handbook of Optical Constant of solid version 1.0 HOC
%
%    1.28 |         | .000124  .04592  | JPH:   Electronic Handbook of Optical Constant of solid version 1.0 HOC
%    1.35 |         | .04592   .3647   | JPH:   Electronic Handbook of Optical Constant of solid version 1.0 HOC
%    1.33 |         | .3757    2.066   | JPH:   Electronic Handbook of Optical Constant of solid version 1.0 HOC
%    1.14 |         | 1.265    9.919   | JPH:   Electronic Handbook of Optical Constant of solid version 1.0 HOC
%    1.5  |         |  .75     .9      | modele de Drude benchmark
%    1.7  |         |  0.1879,1.9372   | P.B. Johnson. R.W. Christy. Phys. Rev. B6. 4370 (1972)
%    1.71 |         |  (        )      | modele de Drude donné par P.B. Johnson. R.W. Christy transmis par Christophe Sauvan
%    1.72 |         |  0.1879,2.0664   | Modele de Lorentz à 4 termes Normander 'Effective dielectric function for FDTD...' Chemical Physics Letters 446 (2007) 115-118 Transmis par JJG
%    1.8  |         |  (        )      | modele de Drude donné par A Archambault,F Marquier, JJ Greffet PRB 82 ,035411 (2010)
%    1.91 |         |  (.25,12.4)      | Modele de Lorentz-Drude A.D.Rakic,A B.Djurisic,J.M.Elazar,M.L.Majewski Appl.Opt.V37,5271-5283(1998)
%    1.92 |         |  (.25,12.4)      | Modele de Brendel_Bormann A.D.Rakic,A B.Djurisic,J.M.Elazar,M.L.Majewski Appl.Opt.V37,5271-5283(1998) 
%    1.61 |         |                  | modele de Drude Ordal transmis par P Ben-Abdallah  correction JPH 2 2014 d'apres le  Palik pour grand ld
%--------------------------------------------------------------------------------------------------------------------------
%    2    |  or     | .15      15      | P.Lalanne:  Palik et Electronic Handbook of Optical Constant of solid version 1.0 HOC
%
%    2.01 |         | .5166    .9537   | JT:    M.C. Theye, Phys. Rev. B2, 3060 (1970)
%    2.03 |         | .5       .95     | JT:    K. Weiss. Z. Naturforscher 3a. 143 (1948)	
%    2.04 |         | .55      .95     | JT:    Corn   http://www.sopra   sa.com/index2.htm 
%    2.05 |         | .51662   .99191  | JT:    http://www.sopra   sa.com/index2.htm
%
%    2.28 |         | .000125  .0124   | JPH:   Electronic Handbook of Optical Constant of solid version 1.0 HOC
%    2.39 |         | .01409   .04768  | JPH:   Electronic Handbook of Optical Constant of solid version 1.0 HOC
%    2.10 |         | .047     .2      | JPH:   Electronic Handbook of Optical Constant of solid version 1.0 HOC
%    2.27 |         | .2066    2.066   | JPH:   Electronic Handbook of Optical Constant of solid version 1.0 HOC
%    2.14 |         | 1.265    4       | JPH:   Electronic Handbook of Optical Constant of solid version 1.0 HOC
%
%    2.2  |         | .20664   2.479   | H.Sauer:    ref='???????? ';
%    2.3  |         | .5       1       |  approximation polynomiale de degre 4 de Corn   (2.04)
% 
%    2.5  |         |  (.5 1.1    ^)   | modele de Drude Lorentz: F.Kaminski,V.Sandoghdar,and M.Agio Journal of Computional and Theoretical Nanoscience Vol 4 635-643 2007
%
%    2.7  |         |  0.1879,1.9372   | P.B. Johnson. R.W. Christy. Phys. Rev. B6. 4370 (1972)
%    2.72 |         |  0.1879,2.0664   | Modele de Lorentz à 4 termes Normander 'Effective dielectric function for FDTD...' Chemical Physics Letters 446 (2007) 115-118 Transmis par JJG
%    2.91 |         |  (.25, 6.2)      | Modele de Lorentz-Drude A.D.Rakic,A B.Djurisic,J.M.Elazar,M.L.Majewski Appl.Opt.V37,5271-5283(1998)
%    2.92 |         |  (.25, 6.2)      | Modele de Brendel_Bormann A.D.Rakic,A B.Djurisic,J.M.Elazar,M.L.Majewski Appl.Opt.V37,5271-5283(1998) 
%    2.61 |         |                  |   correction JPH 2 2014 d'apres le  Palik pour grand ld
%--------------------------------------------------------------------------------------------------------------------------
%    3    | Chrome  | 0.20664   1.239  | 
%   3.07  | Chrome  |  0.1879,1.9372   | P.B. Johnson. R.W. Christy. Phys. Rev. B6. 5059
%    3.26 | Chrome  | 0.01   20 |Palick  et Electronic Handbook of Optical Constant of solid version 1.0 HOC transmis par Mathieu
%    3.7  |         | 0.0100     30.995| ????  Concatenation of: Cr_llnl_cxro + Cr_palik   http://www.luxpop.com????? pas compatible avec  3 et 20.7
%    3.91 |         |  (.25,  62)      | Modele de Lorentz-Drude A.D.Rakic,A B.Djurisic,J.M.Elazar,M.L.Majewski Appl.Opt.V37,5271-5283(1998)
%    3.92 |         |  (.25,  62)      | Modele de Brendel_Bormann A.D.Rakic,A B.Djurisic,J.M.Elazar,M.L.Majewski Appl.Opt.V37,5271-5283(1998) 
%--------------------------------------------------------------------------------------------------------------------------
%    4.1  | SF11    |                  | 
%    4.2  | SF10    |                  | 
%    4.3  | BK7     |                  |
%    4.31 | SiO2    |                  | Sellmeier equation Wikipedia ajout 1 2013
%    4.4  | SiO2    | 3.7070, .2130    | extrait de la bibliothèque CodeV
%    4.41 | SiO2    | 1.5385, 20       | Palik transmis par Simon
%    4.42 | SiO2    |                  | modele de Lorentz d'apres 4.41 fité par JPH valable dans le domaine  1.5385, 20
%    4.45 |MoSiO2   |0.0001239   .41326|  Sébastien de Rossi. Calculé à partir des 'atomic scattering factors' de (6.31g/cm3) CXRO <http://www-cxro.lbl.gov> et LLNL <http://www-phys.llnl.gov/V_Div/scattering/asf.html>.
%    4.5  |Silicium |   1.1    5.6     | Frey Leviton Madison,Proc,SPIE Vol,6273,6272J'  
%    4.53 |Silicium |   .25    1.45    | Green,M.A.,and M.J.Keevers,Optical properties of intrinsic silicon at 300 K,Progress in Photovoltaics:Research and Applications,vol.3,issue 3,pp.189-192,1995 (transmis par Anthony Jouanin)
%    4.52 |Silicium |0.0001239   3332. |  Sébastien de Rossi      Si_llnl_cxro + Si_palik   
%    4.51 |Germanium|   1.9    5.5     | Frey Leviton Madison,Proc,SPIE Vol,6273,6272J'  
%    4.6  |AlxGA1-xAs|                 | Optoélectronique, Emmanuel Rosencher et Borge Vinter, Masson, chapitre 12.7, P332'; 
%    4.61 |PxAs1-xGa|                  | Optoélectronique, Emmanuel Rosencher et Borge Vinter, Masson, chapitre 12.7, P332'; 
%    4.62 |GaxIn1-xP|                  | Optoélectronique, Emmanuel Rosencher et Borge Vinter, Masson, chapitre 12.7, P332'; 
%    4.63 |GaAs dopé|                  | formule de  Simon Vassant Palik  Electronic Handbook of Optical Constant of solid ISBN 0-12-544420-60 p429
%    4.71 |Saphir no| 1.1098 1.6932    | Laurent FREY Handbook of Optical   Constants of Solids III, Ed. D. Palik (1998) transmis par P Lalanne
%    4.72 |Saphir ne|                  |
%    4.8  |TiO2     |.35            .8 | 'Yamada Deposition at low substrate temperatures of hight quality TiO2 Films by radical beam assistef evaporation Appl Opt 1999 38 6638-6641 courbe Fig2 300K fitee par JPH ';
%   4.81  |TiO2     | 0.001      4.5   |Palick  et Electronic Handbook of Optical Constant of solid version 1.0 HOC, transmis par mathieu
%   4.82  |TiO2     |                  |Sellmeier equation transmis par  Mathieu http://refractiveindex.info
%   4.9   |MgF2     |                  |Sellmeier equation transmis par Mathieu  http://refractiveindex.info   
%  4.11   |Al2O3    |                  |Sellmeier equation transmis par Mathieu  http://refractiveindex.info   
%  4.14   |SiO2     |                  |Sellmeier equation transmis par Mathieu  http://refractiveindex.info
%  4.12   |ZnO      |                  |Sellmeier equation transmis par Mathieu  http://refractiveindex.info
%  4.15   | Si3N4   |0.18    4.8       | Palik (?) transmis par Mathieu';
%--------------------------------------------------------------------------------------------------------------------------
%    5    | eau     |                  | 
%--------------------------------------------------------------------------------------------------------------------------
%    10   |aluminium| 0.00415  30.99   | http://cvs.savannah.nongnu.org/viewvc/freesnell/al.nk?root=freesnell&view=markup  project freesnell - View of /freesnell/al.nk 23 Octobre 2004 (Transmis par Inbal Friedler et P.Lalanne)
%    10.1 |         |                  | TRES MAUVAIS *** modele de Drude PhysRevB.68.245405
%   10.91 |         |  (        )      | Modele de Lorentz-Drude A.D.Rakic,A B.Djurisic,J.M.Elazar,M.L.Majewski Appl.Opt.V37,5271-5283(1998)
%   10.92 |         |  (        )      | Modele de Brendel_Bormann A.D.Rakic,A B.Djurisic,J.M.Elazar,M.L.Majewski Appl.Opt.V37,5271-5283(1998) 
%--------------------------------------------------------------------------------------------------------------------------
%    12   |  cuivre |.000137756,9.53692| |  A VERIFIER | Handbook of Optical Constants of Solids, Ed. by Edward D. Palik,http://www.luxpop.com
%    12.7 |         |  0.1879,1.9372   | P.B. Johnson. R.W. Christy. Phys. Rev. B6. 4370 (1972)
%   12.91 |         |  (        )      | Modele de Lorentz-Drude A.D.Rakic,A B.Djurisic,J.M.Elazar,M.L.Majewski Appl.Opt.V37,5271-5283(1998)
%   12.92 |         |  (        )      | Modele de Brendel_Bormann A.D.Rakic,A B.Djurisic,J.M.Elazar,M.L.Majewski Appl.Opt.V37,5271-5283(1998) 
%   12.61 |         |                  | modele de Drude Ordal transmis par P Ben-Abdallah  correction JPH 2 2014
%--------------------------------------------------------------------------------------------------------------------------
%    13   | Tungsten|0.000012398,24.7960|  'W_llnl_cxro + W_palik'; http://www.luxpop.com (corrigé en 11 2013 )
%  13.24  |Tungsten | 0.1       24.80  | Electronic Handbook of Optical Constant of solid version 1.0 HOC transmis par Mathieu
%  13.241 |         | 0.1       4.5    |Electronic Handbook of Optical Constant of solid version 1.0 HOC transmis par Mathieu
%   13.91 |         |  (.25, 12.4)     | Modele de Lorentz-Drude A.D.Rakic,A B.Djurisic,J.M.Elazar,M.L.Majewski Appl.Opt.V37,5271-5283(1998)
%   13.92 |         |  (.25, 12.4)     | Modele de Brendel_Bormann A.D.Rakic,A B.Djurisic,J.M.Elazar,M.L.Majewski Appl.Opt.V37,5271-5283(1998) 
%--------------------------------------------------------------------------------------------------------------------------
%    14.7 | Fer     |  0.1879,1.9372   | P.B. Johnson. R.W. Christy. Phys. Rev. B6. 5059
%--------------------------------------------------------------------------------------------------------------------------
%    15.7 | Cobalt  |  0.1879,1.9372   | P.B. Johnson. R.W. Christy. Phys. Rev. B6. 5059
%--------------------------------------------------------------------------------------------------------------------------
%    16.7 | Nickel  |  0.1879,1.9372   | P.B. Johnson. R.W. Christy. Phys. Rev. B6. 5059
%   16.91 |         |  (.25, 6.2)      | Modele de Lorentz-Drude A.D.Rakic,A B.Djurisic,J.M.Elazar,M.L.Majewski Appl.Opt.V37,5271-5283(1998)
%   16.92 |         |  (.25, 6.2)      | Modele de Brendel_Bormann A.D.Rakic,A B.Djurisic,J.M.Elazar,M.L.Majewski Appl.Opt.V37,5271-5283(1998) 
%--------------------------------------------------------------------------------------------------------------------------
%    17.7 |Palladium|  0.1879,1.9372   | P.B. Johnson. R.W. Christy. Phys. Rev. B6. 5059
%   17.91 |         |  (.25,12.4)      | Modele de Lorentz-Drude A.D.Rakic,A B.Djurisic,J.M.Elazar,M.L.Majewski Appl.Opt.V37,5271-5283(1998)
%   17.92 |         |  (.25,12.4)      | Modele de Brendel_Bormann A.D.Rakic,A B.Djurisic,J.M.Elazar,M.L.Majewski Appl.Opt.V37,5271-5283(1998) 
%--------------------------------------------------------------------------------------------------------------------------
%    18   |  Titane | 0.000123  200    | Electronic HandBook of Optical Constants of Solids (HOC) transmis par JC Rodier
%    18.7 |         |  0.1879,1.9372   | P.B. Johnson. R.W. Christy. Phys. Rev. B6. 5059
%   18.91 |         |  (.25,  30)      | Modele de Lorentz-Drude A.D.Rakic,A B.Djurisic,J.M.Elazar,M.L.Majewski Appl.Opt.V37,5271-5283(1998)
%   18.92 |         |  (.25,  30)      | Modele de Brendel_Bormann A.D.Rakic,A B.Djurisic,J.M.Elazar,M.L.Majewski Appl.Opt.V37,5271-5283(1998) 
%   18.61 |         |                  |   modele de Drude Ordal transmis par P Ben-Abdallah  correction JPH 2 2014
%--------------------------------------------------------------------------------------------------------------------------
%    19.7 | Vanadium|  0.1879,1.9372   | P.B. Johnson. R.W. Christy. Phys. Rev. B6. 5059
%--------------------------------------------------------------------------------------------------------------------------
%    21.7 |Manganese|  0.1879,1.9372   | P.B. Johnson. R.W. Christy. Phys. Rev. B6. 5059
%--------------------------------------------------------------------------------------------------------------------------
%    22   |MgO      |  12.2, 43.4      | Handbook of Optical Constants of Solids, Ed. by Edward D. Palik, transmis par Mondher
%    22.1 |MgO      | 0.01653   20     | Palick  et Electronic Handbook of Optical Constant of solid
%--------------------------------------------------------------------------------------------------------------------------
%    23   |SiC      |                  | modele de Lorentz donnees Palik transmis par:J le Gall M Olivier JJ Greffet Phys Rev B vol 55 number 15 15 april 1997
%    23.1 |SiC      |0.04133,   25     | Palik hexagonal SiC, ordinary ray.  L'accord avec 23 est remarquable dans le domaine 5 25 mais Palik manque de points dans la zone 12.4 12.8: 23 est preferable 
%--------------------------------------------------------------------------------------------------------------------------
%   24    |Platine  | 0.3757  12.4     | Palik transmis par Daniele Costantini
%   24.91 |         |  (.25, 12.4)     | Modele de Lorentz-Drude A.D.Rakic,A B.Djurisic,J.M.Elazar,M.L.Majewski Appl.Opt.V37,5271-5283(1998)
%   24.92 |         |  (.25, 12.4)     | Modele de Brendel_Bormann A.D.Rakic,A B.Djurisic,J.M.Elazar,M.L.Majewski Appl.Opt.V37,5271-5283(1998) 
%--------------------------------------------------------------------------------------------------------------------------
%   25.91 |Beryllium|  (.25,  62)      | Modele de Lorentz-Drude A.D.Rakic,A B.Djurisic,J.M.Elazar,M.L.Majewski Appl.Opt.V37,5271-5283(1998)
%   25.92 |         |  (.25,  62)      | Modele de Brendel_Bormann A.D.Rakic,A B.Djurisic,J.M.Elazar,M.L.Majewski Appl.Opt.V37,5271-5283(1998)
%--------------------------------------------------------------------------------------------------------------------------
%    27   | Tantale | 0.01      125    | Palik, transmis par Mathieu
%--------------------------------------------------------------------------------------------------------------------------
%    26   |Molybdene| 0.01      12.40  |Electronic Handbook of Optical Constant of solid version 1.0 HOC, transmis par Mathieu
%    26.1 |         | 0.01      4.5    |Electronic Handbook of Optical Constant of solid version 1.0 HOC, transmis par Mathieu
%-------------------------------------------------------------------------
% 
%  method :methode d'interpolation dans le cas de tables. 
% ATTENTION Par defaut:method='spline' (voir interp1)
% il est plus prudent d'utiliser 'linear' qui conserve le signe de la partie imaginaire !
% LES VERSIONS 'OFFICIELLES' DE L'Or ET L'ARGENT SONT 2 ET 1 avec method='spline' (option par defaut) 
%
%  dans le cas 4.5 ou 4.51 (Silicium Germanium) peut etre remplacee par la temperature en Kelvin  (293 par defaut)
%  dans le cas 4.6 (AlxGA1-xAs) peut etre remplacee par x  (0 par defaut)
% dans le cas 4.63 peut etre remplacee par Nd  (0 par defaut)  ou [Nd,Gamma_e] Gamma_e=60 par defaut
% 	retindice(35,4.63) GaAs pas dopé
% 	retindice(35,4.63,1.5e6) GaAs dopé Nd=1.5e6 /mu^3=1.5e18 /cm^3 Gamma_e=60
% 	retindice(35,4.63,[1.5e6,40]) GaAs dopé Nd=1.5e6 /mu^3=1.5e18 /cm^3 Gamma_e=40
% 
%% TESTS DES DONNEES: 
%  si pas de variable en sortie, trace de la courbe sauf si ld a une seule valeur
%  si parm est un tableau, trace comparatif des paramatres
%     si en plus ll=[],trace de la courbe entre les bornes les plus grandes des parametres choisis 
%   (dans le cas de tabulation,les points tabules sont marqués) 
%    si parm a une partie imaginaire, trace de epsilon=n.^2
% Remarque Les tableaux sont dans le fichier retindice.mat et sont lus et mis dans des variables persistent à la premiére utilisation
%
% formes particulieres:
%----------------------------
% n=retindice(ld,d,n2,n1) :indice effectif du mode fondamental transverse d'une fibre circulaire de diamétre d d'indice n2 dans un milieu d'indice n1
% d et ld ont la meme longueur,(ou bien l'un des 2 a un seul element)  n2 et n1 sont reels   n2/n1 peut varier entre 1.18 et 5
%   ( principe: on part d'une formule approchée du type: (n-n1)/(n2-n1) = (1-1/(a3*x^2+a4)) * exp(-exp(a1*x+a2)) , x=n2*d./ld , 
%     suivie d'une recherche précise avec retmarcuse. retindice retourne bien sur ce dernier résultat)
%
%  n=retindice(ld,ld0,ctho,n_ext,h_boite,r_boite,F);
%  n_ext indice exterieur
%  ld0,ctho longueur d'onde de resonance, distance parcourue par la lumiere dans le vide pendant la duree de vie (microns)
%  h_boite,r_boite: hauteur et rayon de la boite (microns)
%  F force de l'oscillateur (facultatif. Par defaut 3F=1)
%
%% EXEMPLES de tracés:
%  retindice([],1);
%  retindice(linspace(.2,2,1000),[2.2,2.27,2.14,2,2.4]);title('Or');
%  retindice([],[2.06,2.2,2.27,2.14,2,2.01,2.02,2.03,2.04,2.05,2.7]);title('Or');
%  retindice(linspace(.5,1,1000),[2.2,2.27,2.14,2,2.01,2.02,2.03,2.04,2.05,2.3]);title('Or');
%  retindice([],[1,1.28,1.35,1.33,1.14,1.7]);title('ARGENT');
%  retindice([],10);title('ALUMINIUM');
%  retindice([],3);
%
% See also: RETHELP_POPOV RETDRUDE RETINDICE_INIT
%
% Internet:   http://refractiveindex.info      http://www.luxpop.com 

if nargin<3;method='spline';end;if nargin<2;parm=1;end;
if nargin==4;[n,bornes,ref]=indicefibre(ld,parm,method,varargin{:});[ldr,nr,ldi,ni]=deal([]);return;end;
if nargin>=6;[n,bornes,ref]=indice_BQ(ld,parm,method,varargin{:});[ldr,nr,ldi,ni]=deal([]);return;end;


%  COMPARAISONS avec superposition de courbes
trace_eps=~all(imag(parm)==0);parm=real(parm);
if length(parm)>1; 
m=length(parm);n=cell(1,m);ldr=cell(1,m);ldi=cell(1,m);nr=cell(1,m);ni=cell(1,m);
if isempty(ld);ld=[inf,-inf];% construction de ld
for ii=1:m;[prv,bornes]=retindice([],parm(ii));if ~isfinite(bornes(2));bornes(2)=10;end;if bornes(1)<=0;bornes(1)=.2;end;ld(1)=min(ld(1),bornes(1));ld(2)=max(ld(2),bornes(2));end;
ld=linspace(ld(1),ld(2),10000);
end;
reel=1;for ii=1:m;[n{ii},bornes,ref,ldr{ii},nr{ii},ldi{ii},ni{ii}]=retindice(ld,parm(ii),method);if ~isreal(n{ii});reel=0;end;end;
ma={'.k','*r','xb','+g','oc','*m','sm','dc','vg','^b','>r','<k'};maa={'-k','--r','-b','--g','c','m','-m','--c','-g','--b','r','k'};
figure;hold on;leg='lg=legend(';
for ii=1:m;jj=1+mod(ii-1,min(length(ma),length(maa)));
f=find(isfinite(n{ii})&~isnan(n{ii}));n{ii}=n{ii}(f);LD=ld(f);
f=find((ldr{ii}<max(ld))&(ldr{ii}>min(ld)));nr{ii}=nr{ii}(f);ldr{ii}=ldr{ii}(f);
f=find((ldi{ii}<max(ld))&(ldi{ii}>min(ld)));ni{ii}=ni{ii}(f);ldi{ii}=ldi{ii}(f);

if trace_eps;
    if reel==1;
    plot(LD,n{ii},maa{jj},ldr{ii},nr{ii},maa{jj});grid on;ylabel('epsilon');    
    else;    
    subplot(2,1,1);hold on;plot(LD,real(n{ii}.^2),maa{jj},'linewidth',2);grid on;ylabel('real(epsilon)');
    subplot(2,1,2);hold on;plot(LD,imag(n{ii}.^2),maa{jj},'linewidth',2);grid on;xlabel('longueur d''onde en microns');ylabel('imag(epsilon)');
    end;
else    
    if reel==1;
    plot(LD,n{ii},maa{jj},ldr{ii},nr{ii},ma{jj});grid on;;ylabel('n');    
    else;    
    subplot(2,1,1);hold on;plot(LD,real(n{ii}),maa{jj},ldr{ii},nr{ii},ma{jj},'linewidth',2);grid on;ylabel('real(n)');
    subplot(2,1,2);hold on;plot(LD,imag(n{ii}),maa{jj},ldi{ii},ni{ii},ma{jj},'linewidth',2);grid on;xlabel('longueur d''onde en microns');ylabel('imag(n)');
    end;
end;

if ~isempty(LD);leg=[leg,'''',num2str(parm(ii)),''','];end;
if ~isempty(ldr{ii}) & ~trace_eps;leg=[leg,''''','];end;
end;

leg=[leg(1:end-1),');'];eval(leg);set(lg,'fontsize',7);xlabel('longueur d''onde en microns');
retfont(gcf,0);

return
end;

% CALCUL
switch parm;
case {1,2};   [n,bornes,ref,ldr,nr,ldi,ni]=indice_philippe(ld,parm,method);  
case {1.14,1.28,1.33,1.35}; [n,bornes,ref,ldr,nr,ldi,ni]=indice_jp(ld,parm,method);  
case {2.10,2.14,2.27,2.28,2.39}; [n,bornes,ref,ldr,nr,ldi,ni]=indice_jp(ld,parm,method); 
case {2.01,2.02,2.03,2.04,2.05}; [n,bornes,ref,ldr,nr,ldi,ni]=indice_jt(ld,parm,method);
case {1.91,1.92,2.91,2.92,12.91,12.92,10.91,10.92,25.91,25.92,3.91,3.92,16.91,16.92,17.91,17.92,24.91,24.92,18.91,18.92,13.91,13.92};[n,bornes,ref,ldr,nr,ldi,ni]=indice_Rakic_metaux(ld,parm);
case {1.6,10.6,2.6,15.6,12.6,11.6,18.6,2.61,1.61,12.61,18.61};[n,bornes,ref,ldr,nr,ldi,ni]=indice_drudePBA(ld,parm);
case 4.7;  [n,bornes,ref,ldr,nr,ldi,ni]=indice_GaAs(ld,parm,method);
case {1.7,2.7,12.7,14.7,15.7,16.7,17.7,18.7,19.7,20.7,3.07,21.7};[n,bornes,ref,ldr,nr,ldi,ni]=indice_Johnson(ld,parm,method);
case {1.72,2.72};[n,bornes,ref,ldr,nr,ldi,ni]=indice_Nordlander(ld,parm);
case 1.71; [n,bornes,ref,ldr,nr,ldi,ni]=indiceag_drude_Johnson(ld);
case 2.06; [n,bornes,ref,ldr,nr,ldi,ni]=indice_g(ld,parm,method); 
case 2.2; [n,bornes,ref,ldr,nr,ldi,ni]=calcul_indice('or',ld,method);  
case 2.3; [n,bornes,ref,ldr,nr,ldi,ni]=indiceor(ld);
case 2.4; [n,bornes,ref,ldr,nr,ldi,ni]=indiceor_drude(ld);
case 2.5; [n,bornes,ref,ldr,nr,ldi,ni]=indiceor_drude_lorentz(ld);
case 1.4; [n,bornes,ref,ldr,nr,ldi,ni]=indiceag_drude(ld);
case 1.5; [n,bornes,ref,ldr,nr,ldi,ni]=indiceag_drude_benchmark(ld);
case 1.8; [n,bornes,ref,ldr,nr,ldi,ni]=indiceag_drude_Alex(ld);
case 4;   [n,bornes,ref,ldr,nr,ldi,ni]=calcul_indice('verre',ld,method);    
case 3;   [n,bornes,ref,ldr,nr,ldi,ni]=calcul_indice('chrome',ld,method); 
case 3.7;[n,bornes,ref,ldr,nr,ldi,ni]=indice_Cr(ld,parm,method);  
case 4;   [n,bornes,ref,ldr,nr,ldi,ni]=calcul_indice('verre',ld,method);    
case 4.1; [n,bornes,ref,ldr,nr,ldi,ni]=calcul_indice('SF11',ld,method);    
case 4.2; [n,bornes,ref,ldr,nr,ldi,ni]=calcul_indice('SF10',ld,method);    
case 4.3; [n,bornes,ref,ldr,nr,ldi,ni]=calcul_indice('BK7',ld,method);    
case 4.31; [n,bornes,ref,ldr,nr,ldi,ni]=calcul_indice('Silice',ld,method);    
case 4.4; [n,bornes,ref,ldr,nr,ldi,ni]=calcul_indice('SiO2',ld,method);   
case 4.41; [n,bornes,ref,ldr,nr,ldi,ni]=calcul_indice('SiO2_Simon',ld,method);
case 4.42;[n,bornes,ref,ldr,nr,ldi,ni]=indiceSiO2_Lorentz(ld);    
case 4.45;[n,bornes,ref,ldr,nr,ldi,ni]=calcul_indice('MoSiO2',ld,method);   
case 4.5; [n,bornes,ref,ldr,nr,ldi,ni]=indice_Si(ld,method,0);   
case 4.52;[n,bornes,ref,ldr,nr,ldi,ni]=calcul_indice('Si',ld,method);  
case 4.53;[n,bornes,ref,ldr,nr,ldi,ni]=calcul_indice('Silicium',ld,method);  
case 4.51;[n,bornes,ref,ldr,nr,ldi,ni]=indice_Si(ld,method,1);   
case 4.6; [n,bornes,ref,ldr,nr,ldi,ni]=indice_algaas(0,ld,method);   
case 4.61;[n,bornes,ref,ldr,nr,ldi,ni]=indice_algaas(1,ld,method);   
case 4.62;[n,bornes,ref,ldr,nr,ldi,ni]=indice_algaas(2,ld,method);   
case 4.63;[n,bornes,ref,ldr,nr,ldi,ni]=indice_GaAs_Simon(ld,method);   
case 4.71;[n,bornes,ref,ldr,nr,ldi,ni]=calcul_indice('Saphir_no',ld,method);   
case 4.8;[n,bornes,ref,ldr,nr,ldi,ni]=calcul_indice('TiO2',ld,method);
case 4.82;[n,bornes,ref,ldr,nr,ldi,ni]=calcul_indice('TiO2_bis',ld,method);
case 4.9;[n,bornes,ref,ldr,nr,ldi,ni]=calcul_indice('MgF2',ld,method);
case 4.11;[n,bornes,ref,ldr,nr,ldi,ni]=calcul_indice('Al2O3',ld,method);
case 4.14;[n,bornes,ref,ldr,nr,ldi,ni]=calcul_indice('SiO2_bis',ld,method);
case 4.12;[n,bornes,ref,ldr,nr,ldi,ni]=calcul_indice('ZnO',ld,method);
case 4.9;[n,bornes,ref,ldr,nr,ldi,ni]=calcul_indice('MgF2',ld,method);   
case 4.11;[n,bornes,ref,ldr,nr,ldi,ni]=calcul_indice('Al2O3',ld,method);   
case 4.12;[n,bornes,ref,ldr,nr,ldi,ni]=calcul_indice('Zn0',ld,method);   
case 4.72;[n,bornes,ref,ldr,nr,ldi,ni]=calcul_indice('Saphir_ne',ld,method);   
case 5;   [n,bornes,ref,ldr,nr,ldi,ni]=calcul_indice('eau',ld,method);    
case 5.1; [n,bornes,ref,ldr,nr,ldi,ni]=calcul_indice('eauSD',ld,method);    
case 6;   [n,bornes,ref,ldr,nr,ldi,ni]=calcul_indice('air',ld,method);    
case 7;   [n,bornes,ref,ldr,nr,ldi,ni]=calcul_indice('bio',ld,method);    
case 8;   [n,bornes,ref,ldr,nr,ldi,ni]=calcul_indice('résine',ld,method);    
case 9;   [n,bornes,ref,ldr,nr,ldi,ni]=calcul_indice('polypyrrole',ld,method);    
case 10;  [n,bornes,ref,ldr,nr,ldi,ni]=calcul_indice('Aluminium',ld,method);
case 10.1; [n,bornes,ref,ldr,nr,ldi,ni]=indicealu_drude(ld);
case {11,18};  [n,bornes,ref,ldr,nr,ldi,ni]=calcul_indice('Titane',ld,method);    
case 12;  [n,bornes,ref,ldr,nr,ldi,ni]=calcul_indice('Cuivre',ld,method);    
case 13;  [n,bornes,ref,ldr,nr,ldi,ni]=calcul_indice('Tungsten',ld,method);    
%case {3.26,4.15,4.81,13.24,13.241,22.1,25.1,26,26.1};  [n,bornes,ref,ldr,nr,ldi,ni]=indice_mathieu(ld,parm,method);    
case 22;  [n,bornes,ref,ldr,nr,ldi,ni]=calcul_indice('MgO',ld,method);    
case 24;  [n,bornes,ref,ldr,nr,ldi,ni]=calcul_indice('Platine',ld,method);    
case 23;  [n,bornes,ref,ldr,nr,ldi,ni]=indiceSiC_drude(ld);
case 23.1;[n,bornes,ref,ldr,nr,ldi,ni]=calcul_indice('SiC',ld,method); 
case {23.2,23.3}; [n,bornes,ref,ldr,nr,ldi,ni]=indiceSi_Lorentz(ld,parm);
case  {28.1,28.2,28.31,28.32,4.13,28.41,28.42,28.43,28.44,4.91,4.701,4.601,4.602,4.6021,4.603,4.604,4.605,4.606};[n,bornes,ref,ldr,nr,ldi,ni]=indice_Collin(ld,parm,method);
case  {13.24,13.241,26,26.1,4.15,4.81,3.26,22.1,27,23.12,23.11,30,31,32,33,40,41,42,43,44};[n,bornes,ref,ldr,nr,ldi,ni]=indice_mathieu(ld,parm,method);
case {101.7,102.7,103.7,104.7,114,114.1,115,116,117,117.1,118,118.1,119,120,114.3,114.2,104.8,115.1,121,118.2,118.3,122};[n,bornes,ref,ldr,nr,ldi,ni]=indice_Collin_interpole(ld,parm,method);
	
otherwise;error(['parm=',num2str(parm),'  inconnu']);
end;

% TRACE 
if (nargout==0)&(length(ld)~=1);
if isempty(ld);if ~isfinite(bornes(2));bornes(2)=15;end;ld=linspace(max(0,bornes(1)),bornes(2),10000);[n,bornes,ref,ldr,nr,ldi,ni]=retindice(ld,parm,method);end;    
figure;f=find(isfinite(n)&~isnan(n));n=n(f);ld=ld(f);
if trace_eps
     if isreal(n);
    plot(ld,n.^2,'-k','linewidth',2);axis tight;grid;title(['parm=',num2str(parm),'  ',ref],'fontsize',7);xlabel(['longueur d''onde en microns   bornes=',num2str(bornes)]);ylabel('epsilon');    
    else    
    subplot(2,1,1);plot(ld,real(n.^2),'-k','linewidth',2);axis tight;grid;;grid;title(['parm=',num2str(parm),'  ',ref],'fontsize',7);ylabel('real(epsilon)');
    subplot(2,1,2);plot(ld,imag(n.^2),'-k','linewidth',2);axis tight;grid;;grid;xlabel(['longueur d''onde en microns   bornes=',num2str(bornes),' method=',method]);ylabel('imag(epsilon)');
    end;
else;
     if isreal(n);
    plot(ld,n,'-k',ldr,nr,'.k','linewidth',2);axis tight;grid;title(['parm=',num2str(parm),'  ',ref],'fontsize',7);xlabel(['longueur d''onde en microns   bornes=',num2str(bornes)]);ylabel('n');    
    else    
    subplot(2,1,1);plot(ld,real(n),'-k',ldr,nr,'.k','linewidth',2);axis tight;grid;;grid;title(['parm=',num2str(parm),'  ',ref],'fontsize',7);ylabel('real(n)');
    subplot(2,1,2);plot(ld,imag(n),'-k',ldi,ni,'.k','linewidth',2);axis tight;grid;;grid;xlabel(['longueur d''onde en microns   bornes=',num2str(bornes),' method=',method]);ylabel('imag(n)');
    end;
end;
retfont(gcf,0);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [n,bornes,ref,ldr,nr,ldi,ni]=indice_philippe(ld,parm,method);
switch parm;
case 1;% indice de l'argent   Electronic Handbook of Optical Constant of solid version 1.0 HOC
lld=[10 9.919 8.856 7.293 5.636 4.133 3.100 2.296 1.59 1.409 1.265 1.127  1.033    0.9537 .8856 .8266 0.7749 0.7293     .6526 .6199 .5904 .5391 .4959 .4592 0.40 0.3875 0.3757 0.3647 0.3542 .3397 0.3306 .3263 .3179 0.31  0.2988 0.2883 0.2755 0.2638 0.248 .2214 0.2066 0.2]; %en µm
nn=[13.30 13.11  10.69 7.461 4.425 2.446 1.387 0.823 .485 0.421 0.375 0.251  0.226    0.198  .163  .1450 0.1430 0.1480     .140  0.131 .1210 .1290 .1300 .144  .173 0.192  0.20   0.186  0.209 0.259 0.371 0.526  0.932 1.323 1.522  1.502  1.441  1.372  1.298 1.208 1.125  1.072]...
+i*[54.10 53.7   49.40 42.50 34.00 25.10 18.80 14.00 9.57 8.370 7.780 7.670  6.990    6.43   5.95  5.500 5.090  4.74        4.15 3.88   3.66 3.250 2.88  2.56  1.95 1.81   1.67   1.61   1.44  1.12  0.813 0.663  0.504 0.647 0.992  1.19   1.31   1.35   1.35  1.3   1.27   1.24 ];
bornes=[  .2     10   ];ref='Ag  Electronic Handbook of Optical Constant of solid version 1.0 HOC adapté par P Lalanne';
% la valeur a 10 a étée ajoutée par extrapolation (retdrude) en 3 2010 (sinon la derniere valeur etait 9.919)
case 2;% indice de l'or   ref Palik  et Electronic Handbook of Optical Constant of solid version 1.0 HOC celle de Philippe 
lld=[0.15 .16 0.17 .18 .2 .2254 .248 .2755 .3024 .3179 .3444 .3757 0.3875 .4 0.4133 0.4275 0.4428 0.4592 0.4769 0.4959 0.5166 0.5391 0.5636 0.6526 0.6888 0.7293 0.7749 0.8266 0.8856 ...
        0.9537 1.512 2 2.48 2.952 3.444 4 4.592 4.959 5.636 6.199 7.085 8 9.184 9.919 15];% en µm
nn=[1.419 1.483 1.519 1.47 1.427 1.452 1.484 1.648 1.812 1.84 1.766 1.696 1.674 1.658 1.636 1.616 1.562 1.426 1.242 .910 .608 .402 .306 .166 .160 .164 .174 .188 .210 ...
        .236 .537 .85 1.205 1.598 2.046 2.6 3.289 3.748 4.611 5.423 6.937 8.5 10.84 12.24 23]...
+i*[1.102 1.106 1.07 1.085 1.215 1.442 1.636 1.852 1.92 1.904 1.846 1.906 1.936 1.956 1.958 1.94 1.904 1.840 1.790 1.840 2.120 2.54 2.88 3.15 3.80 4.35 4.86 5.39 5.88 ...
        6.47 9.58 12.6 15.5 18.3 21.3 24.6 28.2 30.5 34.3 37.5 42.0 46.4 51.6 54.7 75];
bornes=[ .15       15     ];ref='Au   Palik  et Electronic Handbook of Optical Constant of solid version 1.0 HOC adapté par P Lalanne';
end  % switch parm
if isempty(ld);n=[];else;n=interp1(lld,nn,ld,method,NaN*(1+1i));end;
ldr=lld;ldi=lld;nr=real(nn);ni=imag(nn);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [n,bornes,ref,ldr,nr,ldi,ni]=indice_Cr(ld,parm,method);
% probleme ????????????????????
persistent ld_nCr_niCr;if isempty(ld_nCr_niCr);  load retindice ld_nCr_niCr;end;
lld=ld_nCr_niCr(:,1)/10000;nn=ld_nCr_niCr(:,2)+i*ld_nCr_niCr(:,3);
if isempty(ld);n=[];else;n=interp1(lld,nn,ld,method,NaN*(1+1i));end;
bornes=[min(lld),max(lld)];ref='CR Concatenation of: Cr_llnl_cxro + Cr_palik   http://www.luxpop.com';
ldr=lld;ldi=lld;nr=real(nn);ni=imag(nn);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [n,bornes,ref,ldr,nr,ldi,ni]=indice_GaAs(ld,parm,method);
% indice de  GaAs Palik p 434 440
persistent ld_nGaAs;if isempty(ld_nGaAs);  load retindice ld_nGaAs;end;
lld=ld_nGaAs(:,1);nn=ld_nGaAs(:,2);
if isempty(ld);n=[];else;n=interp1(lld,nn,ld,method,NaN*(1+1i));end;
bornes=[min(lld),max(lld)];ref=['GaAs   Palik p 434 440'];
ldr=lld;ldi=lld;nr=real(nn);ni=imag(nn);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [n,bornes,ref,ldr,nr,ldi,ni]=indice_Johnson(ld,parm,method);
% indice du Cuivre Argent Or fer Cobalt Nickel Palladium Titane Vanadium Chrome  Manganese  P.B. Johnson. R.W. Christy. Phys. Rev. B6. 4370  et  5059
persistent V_nCu_nAG_nOr_nFer_nCo_nNi_nPa_nTi_nVa_nCr_nMg;if isempty(V_nCu_nAG_nOr_nFer_nCo_nNi_nPa_nTi_nVa_nCr_nMg);  load retindice V_nCu_nAG_nOr_nFer_nCo_nNi_nPa_nTi_nVa_nCr_nMg;end;
lld=1.239829488174885./V_nCu_nAG_nOr_nFer_nCo_nNi_nPa_nTi_nVa_nCr_nMg(:,1);
switch parm;
case 12.7;nn=V_nCu_nAG_nOr_nFer_nCo_nNi_nPa_nTi_nVa_nCr_nMg(:,2);ref='Cu ';
case 1.7;nn=V_nCu_nAG_nOr_nFer_nCo_nNi_nPa_nTi_nVa_nCr_nMg(:,3);ref='Ag ';
case 2.7;nn=V_nCu_nAG_nOr_nFer_nCo_nNi_nPa_nTi_nVa_nCr_nMg(:,4);ref='Au ';
case 14.7;nn=V_nCu_nAG_nOr_nFer_nCo_nNi_nPa_nTi_nVa_nCr_nMg(:,5);ref='Fe ';
case 15.7;nn=V_nCu_nAG_nOr_nFer_nCo_nNi_nPa_nTi_nVa_nCr_nMg(:,6);ref='Co  ';
case 16.7;nn=V_nCu_nAG_nOr_nFer_nCo_nNi_nPa_nTi_nVa_nCr_nMg(:,7);ref='Ni ';
case 17.7;nn=V_nCu_nAG_nOr_nFer_nCo_nNi_nPa_nTi_nVa_nCr_nMg(:,8);ref='Pd ';
case 18.7;nn=V_nCu_nAG_nOr_nFer_nCo_nNi_nPa_nTi_nVa_nCr_nMg(:,9);ref='Ti  ';
case 19.7;nn=V_nCu_nAG_nOr_nFer_nCo_nNi_nPa_nTi_nVa_nCr_nMg(:,10);ref='V  ';
case {20.7,3.07};nn=V_nCu_nAG_nOr_nFer_nCo_nNi_nPa_nTi_nVa_nCr_nMg(:,11);ref='Cr ';
case 21.7;nn=V_nCu_nAG_nOr_nFer_nCo_nNi_nPa_nTi_nVa_nCr_nMg(:,12);ref='Mn ';
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(ld);n=[];else;n=interp1(1./lld,nn,1./ld,method,NaN*(1+1i));end;% pour compatibilite avec Collin on interpole en Ev
bornes=[min(lld),max(lld)];ref=[ref,' P.B. Johnson. R.W. Christy. Phys. Rev. B6. 4370  et  5059'];
ldr=lld;ldi=lld;nr=real(nn);ni=imag(nn);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [n,bornes,ref,ldr,nr,ldi,ni]=indice_jp(ld,parm,method);
% indice de l'ARGENT
%  14   B.Dold,R.Mecke.Optik 22,435(1965)
%  28   H.J.Hagemann,W.Gudat,C.Kunz.DESY Report Sr-74/7,Hamburg,(1974) (voir aussi 40 ?)
%  33   P.Winsemius,F.F.van Kampen,H.P.Lengkeek,C.G.van Went.J.Phys.F6,158351976)
%  35   G.Leveque,C.G.Olson,D.W.Lynch.Phys.Rev.B 24,4654 (1983)
persistent ldnag;if isempty(ldnag);  load retindice ldnag;end;
% indice de l'Or
%  10   L.R.Canfield,G.Hass,W.R.Hunter.J.Physique 25,124(1964).
%  14   B.Dold,R.Mecke.Optik 22,435(1965)
%  27   M.L.Theye.Phys.Rev.B2,3060(1970)
%  28   P.O.Nilsson.Phys.Kondens.Mater.11,1(1970)
%  39   H.J.Hagemann,W.Gudat,C.Kunz.DESY Report Sr-74/7,Hamburg,(1974) (voir aussi 55 ? )
persistent ldnor;if isempty(ldnor);  load retindice ldnor;end;
switch parm;
case 1.14;f=find(ldnag(:,3)==14);lld=ldnag(f,1);nn=ldnag(f,2);ref='Ag ';
case 1.28;f=find(ldnag(:,3)==28);lld=ldnag(f,1);nn=ldnag(f,2);ref='Ag ';
case 1.33;f=find(ldnag(:,3)==33);lld=ldnag(f,1);nn=ldnag(f,2);ref='Ag ';
case 1.35;f=find(ldnag(:,3)==35);lld=ldnag(f,1);nn=ldnag(f,2);ref='Ag ';
case 2.10;f=find(ldnor(:,3)==10);lld=ldnor(f,1);nn=ldnor(f,2);ref='Au ';
case 2.14;f=find(ldnor(:,3)==14);lld=ldnor(f,1);nn=ldnor(f,2);ref='Au ';
case 2.27;f=find(ldnor(:,3)==27);lld=ldnor(f,1);nn=ldnor(f,2);ref='Au ';
case 2.28;f=find(ldnor(:,3)==28);lld=ldnor(f,1);nn=ldnor(f,2);ref='Au ';
case 2.39;f=find(ldnor(:,3)==39);lld=ldnor(f,1);nn=ldnor(f,2);ref='Au ';
end;

if isempty(ld);n=[];else;n=interp1(lld,nn,ld,method,NaN*(1+1i));end;
bornes=[min(lld),max(lld)];ref=[ref,'  JPH:   Electronic Handbook of Optical Constant of solid version 1.0 HOC'];
ldr=lld;ldi=lld;nr=real(nn);ni=imag(nn);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [n,bornes,ref,ldr,nr,ldi,ni]=indice_g(ld,parm,method);
% indice de l'Or  transmis par Guillaume
persistent ldn;if isempty(ldn);  load retindice ldn;end;
lld=ldn(:,2)/1000;nn=ldn(:,3)+i*ldn(:,4);
if isempty(ld);n=[];else;n=interp1(lld,nn,ld,method,NaN*(1+1i));end;
bornes=[min(lld),max(lld)];ref=['Au  guillaume sat'];
ldr=lld;ldi=lld;nr=real(nn);ni=imag(nn);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ù%%%%%%%%%%%%%%%%%%%
function [n,bornes,ref,ldr,nr,ldi,ni]=indice_jt(ld,parm,method);
% indice de l'Or  transmis par Jean Taboury
% M.C. Theye, Phys. Rev. B2, 3060 (1970)						
persistent ldnr1 ldni1 ldnr2 ldni2 ldnr3 ldni3 ldnr4 ldni4 ldnr44 ldni44 ldnr5 ldni5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch parm;
case 2.01;
if isempty(ldnr1);  load retindice ldnr1 ldni1;end;% M.C. Theye, Phys. Rev. B2, 3060 (1970)						
lld=ldnr1(:,1);nn=ldnr1(:,2)+i*ldni1(:,2);ref='Au   JT:  M.C. Theye, Phys. Rev. B2, 3060 (1970)	';
case 2.02;
if isempty(ldnr2); load retindice ldnr2 ldni2;end;%P.B. Johnson. R.W. Christy. Phys. Rev. B6. 4370 (1972)
lld=ldnr2(:,1);nn=ldnr2(:,2)+i*ldni2(:,2);ref='Au   JT:     P.B. Johnson. R.W. Christy. Phys. Rev. B6. 4370 (1972)	';
case 2.03;
if isempty(ldnr3);  load retindice ldnr3 ldni3;end;%K. Weiss. Z. Naturforscher 3a. 143 (1948)
lld=ldnr3(:,1);nn=ldnr3(:,2)+i*ldni3(:,2);ref='Au  JT:     K. Weiss. Z. Naturforscher 3a. 143 (1948)	';
case 2.04;% a priori ldnr44 ldni44 nn ne servent pas
if isempty(ldnr4);  load retindice ldnr4 ldni4 ldnr44 ldni44;end;%http://corndog.chem.wisc.edu/fresnel/audata.txt					(Corn)				0.004998306
lld=ldnr44(:,1);nn=ldnr44(:,2)+i*ldni44(:,2);ref='Au   JT:     Corn	 http://corndog.chem.wisc.edu/fresnel/audata.txt';
case 2.05;
if isempty(ldnr5);  load retindice ldnr5 ldni5;end;%http://www.sopra   sa.com/index2.htm
lld=ldnr5(:,1);nn=ldnr5(:,2)+i*ldni5(:,2);ref='Au  JT:   http://www.sopra   sa.com/index2.htm';
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if parm==2.04;
if isempty(ld);n=[];else;n=interp1(ldnr4(:,1),ldnr4(:,2),ld,method,NaN)+i*interp1(ldni4(:,1),ldni4(:,2),ld,method,NaN);end
ldr=ldnr4(:,1);ldi=ldni4(:,1);nr=ldnr4(:,2);ni=ldni4(:,2);
else;
if isempty(ld);n=[];else;n=interp1(lld,nn,ld,method,NaN*(1+1i));end
ldr=lld;ldi=lld;nr=real(nn);ni=imag(nn);
end;
bornes=[min(lld),max(lld)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [indice,bornes,ref,ldr,nr,ldi,ni]= calcul_indice(milieu,lambda,method);
%  calcul l' indice en fonction de lambda pour différents matériaux
%
% USAGE: [indice,bornes,ref]= calcul_indice(milieu,lambda)
% milieu: chaîne de caractères donnant le nom du milieu
% lambda: tableau réel donnant les lambda  EN MICRONS où l'indice doit être calculé
% indice : tableau réel ou complexe de même taille que lambda retournant les indices
%
% bornes:domaine de validite (quand interploation retourne nan a l'exterieur)
% ref:references(chaine de characteres)
%  ldr,nr,ldi,ni   couples de points en cas d'interpolation
%
% LISTE DES MILIEUX
% 'air'
% 'bio'
% 'verre' % NOTA: Sans dispersion (Cf. BK7, SF11, etc...)
% 'résine' % NOTA: Sans dispersion 
% 'eauSD' % Sans Dispersion ('eau' pour le programme IDL initial)
% 'eau' % HS!      Coefficients série de Laurent d'après la Glass Database du catalogue "Misc" de OSLO LT 6.1
% 'TiO2';tmp=lambda.^2;bornes=[-inf,inf];ref='TiO2 Sellmeier equation transmis par Mathieu';    
% 'MgF2';tmp=lambda.^2;bornes=[-inf,inf];ref='MgF2 Sellmeier equation transmis par Mathieu';    
% 'Al2O3';tmp=lambda.^2;bornes=[-inf,inf];ref='Al2O3 Sellmeier equation transmis par Mathieu';    
% 'SiO2';tmp=lambda.^2;bornes=[-inf,inf];ref='SiO2 Sellmeier equation transmis par Mathieu';    
% 'ZnO';tmp=lambda.^2;bornes=[-inf,inf];ref='Al2O3 Sellmeier equation transmis par Mathieu';    
% 'polypyrrole'
% 'SF11'
% 'SF10'
% 'BK7'
% 'Silice'
% 'or','Au'
% 'chrome','Cr'
% 'SiO2'
% 'SiO2_bis'
% 'TiO2'
% 'TiO2_bis'
% 'SiO2_Simon'
% 'MoSiO2'
% 'Si'
% 'Saphir_no'
% 'Saphir_ne'
% 'Aluminium'
% 'Titane'
% 'Cuivre'
% 'Tungsten'
% 'MgO'
% 'Platine'
% 'SiC'
% Une erreur est déclenchée si le milieu est inconnu. 
%
% D'après le programme IDL "determ_indice_2.pro" de Yan Pailhas (d'après les données et la
% biblio de Genoptics et de l'équipe d'Yves Levy...)  - mai-juin 2003
%
% d'apres calcul_epsilon © LCFIO - Hervé Sauer - 8 juin 2003

persistent nn_or nn_chrome AAA_Aluminium AAA_Titane AAA_Cuivre AAA_Tungsten AAA_Si AAA_MoSiO2 ldn_MgO ldn_Silicium AAA_SiO2_Simon AAA_Platine AAA_SiC;
indice=ones(size(lambda));lld=[];nn=[];
switch milieu
  
  case 'air';bornes=[-inf,inf];ref='Air sans dispersion';
  case 'bio';indice=1.45*indice;bornes=[-inf,inf];ref='Bio sans dispersion';
  case 'polypyrrole';indice=1.7*indice;bornes=[-inf,inf];ref='Polypyrrole sans dispersion';
  case 'résine';indice=1.5*indice;bornes=[-inf,inf];ref='Resine sans dispersion';
  case 'eauSD';indice=1.333*indice;bornes=[-inf,inf];ref='H2O sans dispersion ';
  case 'eau' % HS!
    bornes=[-inf,inf];ref='H2O avec dispersion, bornes a preciser,série de Laurent d''après la Glass Database du catalogue "Misc" de OSLO LT 6.1';
    CL=[1.76755000 -0.01507500 0.00350000 0.00061850 -2.820E-5 0];
    % ^ Coefficients série de Laurent d'après la Glass Database du catalogue "Misc" de OSLO LT 6.1
    % v Formule d'après l'aide de CodeV 9.20 (coef pour lambda en µm)
    lb_um=lambda;
    indice=retsqrt(CL(1) + CL(2)*lb_um.^2 + CL(3)*lb_um.^-2 + CL(4)*lb_um.^-4 +CL(5)*lb_um.^-6 +  CL(6)*lb_um.^-8,-1);
    % NOTA: nd=1.333037  d'après OSLO LT 6.1
    % ====: Vd=(nd-1)/(nF-nC)=57.07 d'après OSLO LT 6.1
    %       d(He): 587.5618nm; F(H): 486.1327nm; C(H): 656.2725nm   d'après aide mémoire Schott 1986 
    %  ==> test de la série ==> résultat: OK pour nd (1.3330371346...);  OK  pour Vd (57.07271...)
    
  case 'verre'; indice=1.515*indice;bornes=[-inf,inf];ref='Verre sans dispersion '; 
  
	
	      
  case 'TiO2_bis';tmp=lambda.^2;bornes=[-inf,inf];ref='TiO2 Sellmeier equation, transmis par Mathieu';    
    indice=retsqrt(5.913+ 0.2441*tmp./(tmp-0.0803),-1);
  case 'MgF2';tmp=lambda.^2;bornes=[-inf,inf];ref='MgF2 Sellmeier equation, transmis par Mathieu';    
    indice=retsqrt(1+0.48755108*tmp./(tmp-0.04338408*0.04338408)+0.39875031*tmp./(tmp-0.09461442*0.09461442)+2.3120353*tmp./(tmp-23.793604* 23.793604),-1);
 case 'Al2O3';tmp=lambda.^2;bornes=[-inf,inf];ref='Al2O3 Sellmeier equation, transmis par Mathieu';    
    indice=retsqrt(1+1.4313493*tmp./(tmp-0.0726631*0.0726631)+0.65054713*tmp./(tmp-0.1193242*0.1193242)+5.3414021*tmp./(tmp-18.028251* 18.028251),-1);
 case 'SiO2_bis';tmp=lambda.^2;bornes=[-inf,inf];ref='SiO2 Sellmeier equation,handbook of Optics [ alpha Quartz:Ordinary ray(o)] transmis par Mathieu';    
    indice=retsqrt(1+0.663044*tmp./(tmp-0.060*0.060)+0.517852*tmp./(tmp-0.106*0.106)+0.175912*tmp./(tmp-0.119*0.119)+0.565380*tmp./(tmp-8.844*8.844)+1.675299*tmp./(tmp-20.742*20.742),-1);
 case 'ZnO';tmp=lambda.^2;bornes=[-inf,inf];ref='Al2O3 Sellmeier equation, transmis par Mathieu';    
    indice=retsqrt(2.81418 + 0.87968*tmp./(tmp-0.3042*0.3042)-0.00711*tmp,-1);

% .............. ajout Mathieu	
	
	
  case 'SF11';tmp=lambda.^2;bornes=[-inf,inf];ref='SF11 bornes a preciser';    
    indice=retsqrt(1+ 1.73848403*tmp./(tmp-1.36068604e-2)+3.11168974e-1*tmp./(tmp-6.15960463e-2)+1.17490871*tmp./(tmp-1.21922711e+2),-1);
  case 'SF10';tmp=lambda.^2;bornes=[-inf,inf];ref='SF10 bornes a preciser';
    indice=retsqrt(1+1.61625977*tmp./(tmp-1.27534559e-2)+2.59229334e-1*tmp./(tmp-5.81983954e-2)+1.07762317*tmp./(tmp-1.16607680e+2),-1);
  case 'BK7';tmp=lambda.^2;bornes=[-inf,inf];ref='BK7 bornes a preciser Sellmeier equation Wikipedia';
    indice=retsqrt(1+ 1.03961212*tmp./(tmp-0.00600069867)+0.231792344*tmp./(tmp-0.0200179144)+1.01046945*tmp./(tmp-103.560653),-1);
  case 'Silice';tmp=lambda.^2;bornes=[-inf,inf];ref='SiO2  Silice_fondue bornes a preciser Sellmeier equation Wikipedia';% ajout 1 2013
    indice=retsqrt(1+ 0.6961663*tmp./(tmp-0.00467914826)+.4079426*tmp./(tmp-.0135120631)+.8974794*tmp./(tmp-97.9340025),-1);
% ..........fin ajout Mathieu	
 
	
 case {'or','Au'}
    bornes=[0.20664,2.479];ref='Au ???????? ';
    % q=1.6022E-19; h=6.6261E-34;c=299792458; V=(h*c/(1.e-6*q))./lambda) en ev   lambda en microns
    V=1.239829488174885./lambda;  
    X_V=linspace(0.5,6.0,111); % energies hv en eV...
    lld=1.239829488174885./X_V;
if isempty(nn_or);  load retindice nn_or;end;nn=nn_or;
if isempty(V);indice=[];else;indice=interp1(X_V,nn,V,method,NaN*(1+1i));end; 
% retourne NaN+NaNi en cas d'extrapolation!
case {'SiO2'}   %(extrait de la bibliothèque CodeV 9.61)
bornes=[3.7070,.2130];ref='SiO2 avec dispersion, (extrait de la bibliothèque CodeV 9.61)';
%   lld =[     3.7070   3.3188   2.9306   2.5423   2.1541   1.7659   1.3777    .9894    .6012    .2130];
%   nn=[ 1.39935  1.41111  1.42088  1.42901  1.43574  1.44132  1.44603  1.45055  1.45800  1.53519];
lld =[ 3.7070,3.3188,2.9306,2.5423,2.1541,1.7659,1.3777,.9894,.7,.65,.6,.6012,.55,.5,.45,.4,.35,.3,.25,.2130];
nn=[1.39935,1.41111,1.42088,1.42901,1.43574,1.44132,1.44603,1.45055 ,1.455294,1.456538,1.458041,1.45800,1.459915,1.462331,1.465571,1.470122,1.476896,1.487793,1.507439,1.53519];
if isempty(lambda);indice=[];else;indice=interp1(lld,nn,lambda,method,NaN*(1+1i));end; 
% retourne NaN+NaNi en cas d'extrapolation!
case {'Saphir_no'}   %(extrait de la bibliothèque CodeV 9.61)
bornes=[1.01398,1.69320];ref='Saphir_no Handbook of Optical Constants of Solids III, Ed. D. Palik (1998) Laurent FREY CEA-G/LETI/DOPT/STCO transmis par P Lalanne';
lld=[1.01398,1.12866,1.36728,1.39506,1.52952,1.69320];
nn=[1.75547,1.75339,1.74936,1.74888,1.74660,1.74368];  
if isempty(lambda);indice=[];else;indice=interp1(lld,nn,lambda,method,NaN*(1+1i));end; 
% retourne NaN+NaNi en cas d'extrapolation!
case {'Saphir_ne'}   %(extrait de la bibliothèque CodeV 9.61)
bornes=[1.01398,1.69320];ref='Saphir_ne Handbook of Optical Constants of Solids III, Ed. D. Palik (1998) Laurent FREY CEA-G/LETI/DOPT/STCO transmis par P Lalanne';
lld=[1.01398,1.12866,1.36728,1.39506,1.52952,1.69320];
nn=[1.74794,1.74549,1.74148,1.74101,1.73874,1.73584];  
if isempty(lambda);indice=[];else;indice=interp1(lld,nn,lambda,method,NaN*(1+1i));end; 
% retourne NaN+NaNi en cas d'extrapolation!
case {'chrome','Cr'}
% d'après article: ????????
bornes=[0.20664,1.239];ref='Cr  ???????? ';
% q=1.6022E-19; h=6.6261E-34;c=299792458;V=(h*c/(1.e-6*q))./lambda);% lambda en microns
V=1.239829488174885./lambda;  
X_V=linspace(1,6.0,101); % energies hv en eV...
lld=1.239829488174885./X_V;
if isempty(nn_chrome); load retindice nn_chrome;end;nn=nn_chrome;
    if isempty(V);indice=[];else;indice=interp1(X_V,nn,V,method,NaN*(1+1i));end; 
    % retourne NaN+NaNi en cas d'extrapolation!
%------------------------------------------------------

case {'SiO2_Simon'}
bornes=[1.5385,20];ref='SiO2 Palik transmis par Simon';
if isempty(AAA_SiO2_Simon);  load retindice AAA_SiO2_Simon;end;AAA=AAA_SiO2_Simon;
nn=AAA(:,2)+i*AAA(:,3);
cmm1=1.e4./lambda;
lld=1.e4./AAA(:,1);
if isempty(cmm1);indice=[];else;indice=interp1(AAA(:,1),nn,cmm1,method,NaN*(1+1i));end; 
%------------------------------------------------------

case {'Aluminium'}
bornes=[0.00415,30.99];ref='AL Savannah.';
if isempty(AAA_Aluminium);  load retindice AAA_Aluminium;end;AAA=AAA_Aluminium;
nn=AAA(:,2)+i*AAA(:,3);
% R=abs((nn-1)./(nn+1)).^2;
% retcompare(R,AAA(:,4)) % 4.1814e-004  erreur sur R < 1.6e-3  et < 1.e-4 de ld= .3 à 29.28
X_V=AAA(:,1); % energies hv en eV...
V=1.239829488174885./lambda;  
lld=1.239829488174885./X_V;
if isempty(V);indice=[];else;indice=interp1(X_V,nn,V,method,NaN*(1+1i));end; 
%------------------------------------------------------
case {'Titane'}
bornes=[ 0.000123,200];ref='Ti Electronic HandBook of Optical Constants of Solids (HOC)';
if isempty(AAA_Titane) ; load retindice AAA_Titane;end;AAA=AAA_Titane;
lld=AAA(:,1);nn=AAA(:,2)+i*AAA(:,3);
if isempty(lambda);indice=[];else;indice=interp1(lld,nn,lambda,method,NaN*(1+1i));end; 
%------------------------------------------------------
case {'Cuivre'}
bornes=[ .000137756 ,9.53692];ref='Cu Handbook of Optical Constants of Solids, Ed. by Edward D. Palik,Academic Press, Inc., 1985.';
if isempty(AAA_Cuivre);  load retindice AAA_Cuivre;end;AAA=AAA_Cuivre;
lld=AAA(:,1)*1.e-4;nn=AAA(:,2)+i*AAA(:,3);
if isempty(lambda);indice=[];else;indice=interp1(lld,nn,lambda,method,NaN*(1+1i));end; 
%------------------------------------------------------
case {'Tungsten'}
bornes=[0.000012398,24.7960];ref='W W_llnl_cxro + W_palik';
% 4 valeurs supprimées au voisinage de .0413
if isempty(AAA_Tungsten);  load retindice AAA_Tungsten;end;AAA=AAA_Tungsten;
lld=AAA(:,1)*1.e-4;nn=AAA(:,2)+(1i*AAA(:,3));
if isempty(lambda);indice=[];else;indice=interp1(lld,nn,lambda,method,NaN*(1+1i));end; 
%------------------------------------------------------
case {'Si'}
bornes=[0.00012398,3332];ref='Si llnl cxro + Si palik';
if isempty(AAA_Si);  load retindice AAA_Si;end;AAA=AAA_Si;
lld=AAA(:,1)*1.e-4;nn=AAA(:,2)+i*AAA(:,3);
if isempty(lambda);indice=[];else;indice=interp1(lld,nn,lambda,method,NaN*(1+1i));end; 
%------------------------------------------------------
case {'MoSiO2'}
%  Sébastien de Rossi   .Calculé à partir des 'atomic scattering factors' de  
%  MoSiO2'(6.31g/cm3)
bornes=[0.00012398,.41326];ref='MoSiO2 Sébastien de Rossi  CXRO <http://www-cxro.lbl.gov> et LLNL <http://www-phys.llnl.gov/V_Div/scattering/asf.html>.';
if isempty(AAA_MoSiO2);  load retindice AAA_MoSiO2;end;AAA=AAA_MoSiO2;
lld=AAA(:,1)/10000;nn=AAA(:,2)+i*AAA(:,3);
if isempty(lambda);indice=[];else;indice=interp1(lld,nn,lambda,method,NaN*(1+1i));end; 
%------------------------------------------------------
case {'MgO'}
bornes=[12.2,43.47];ref='MgO Palik transmis par Mondher';
if isempty(ldn_MgO) ; load retindice ldn_MgO;end;ldn=ldn_MgO;
lld=ldn(:,1);nn=ldn(:,2);
if isempty(lambda);indice=[];else;indice=interp1(lld,nn,lambda,method,NaN*(1+1i));end; 
%------------------------------------------------------
case {'Platine'}
bornes=[0.3757,43.47];ref='Pt  Palik transmis par Daniele';
if isempty(AAA_Platine);load retindice AAA_Platine;end;AAA=AAA_Platine;
lld=AAA(:,1);nn=AAA(:,2)+i*AAA(:,3);
if isempty(lambda);indice=[];else;indice=interp1(lld,nn,lambda,method,NaN*(1+1i));end; 
%------------------------------------------------------
case {'SiC'}
bornes=[0.04133,25];ref='SiC Palik hexagonal SiC, ordinary ray http://luxpop.com/Material/SiC_palik.nk';
if isempty(AAA_SiC);load retindice AAA_SiC;end;AAA=AAA_SiC;
lld=AAA(:,1);nn=AAA(:,2)+i*AAA(:,3);
if isempty(lambda);indice=[];else;indice=interp1(lld,nn,lambda,method,NaN*(1+1i));end; 
%------------------------------------------------------
case {'Silicium'}
bornes=[.25,1.45];ref='Si  Green,M.A.,and M.J.Keevers,Optical properties of intrinsic silicon at 300 K,Progress in Photovoltaics:Research and Applications,vol.3,issue 3,pp.189-192,1995 (transmis par Anthony Jouanin)';
%ld=a(:,1)*1.e-3;n=a(:,3)+(1.e-7i/(4*pi))*(a(:,1).*a(:,2));
if isempty(ldn_Silicium);  load retindice ldn_Silicium;end;ldn=ldn_Silicium;
lld=ldn(:,1);nn=ldn(:,2);
if isempty(lambda);indice=[];else;indice=interp1(lld,nn,lambda,method,NaN*(1+1i));end; 
%---------------------------------------------------
case {'TiO2'}   %(extrait de la bibliothèque CodeV 9.61)
bornes=[.35 ,.8];ref='TiO2 Yamada Deposition at low substrate temperatures of hight quality TiO2 Films by radical beam assistef evaporation Appl Opt 1999 38 6638-6641 courbe Fig2 300K fitee par JPH ';
lld =[.35,.4,.45,.5,.6,.7,.8];
nn=[2.92+1.53e-3i,2.6+1.215e-3i,2.49+.99e-3i,2.44+8.19e-4i,2.39+6.3e-4i,2.36+4.95e-4i,2.35+4.05e-4i];
if isempty(lambda);indice=[];else;indice=interp1(lld,nn,lambda,method,NaN*(1+1i));end; 

  otherwise
    error(['Milieu "' milieu '" inconnu!!'])
end
ldr=lld;ldi=lld;nr=real(nn);ni=imag(nn);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [indice,bornes,ref,ldr,nr,ldi,ni]=indice_Collin(lambda,parm,method);
% Transmis par Stehpane Collin 2014

% {28.1,28.2,28.31,28.32,28.41,28.42,23.43,28.44,4.9,4.701,4.601,4.602,4.6021,4.603,4.604,4.605,4.606,4.11,}
persistent CCC_ZnS_cubic  CCC_Si3N4  CCC_MgF2_ordinary  CCC_MgF2_extraordinary  CCC_GaSb ...
    CCC_TiO2_ordinary  CCC_TiO2_extraordinary  CCC_TiO2_SOPRA_1  CCC_TiO2_SOPRA_2  CCC_GaAs_dope_Si ...
    CCC_GaAs  CCC_Al20Ga80As  CCC_Al315Ga685As  CCC_Al30Ga70As  CCC_Al50Ga50As  CCC_Al70Ga30As  CCC_Al80Ga20As  CCC_AlAs;
indice=ones(size(lambda));lld=[];nn=[];

switch parm;
%------------------------------------------------------
case 28.1;% indice de ZnS cubic from Palik
bornes=[ 0.000124 ,   666.7];ref='ZnS_cubic Palik Optical constants I';
if isempty(CCC_ZnS_cubic) ; load retindice CCC_ZnS_cubic;end;AAA=CCC_ZnS_cubic;
%------------------------------------------------------
case 28.2;% indice de Si3N4 from Palik
bornes=[0.05166, 1.24];ref='Si3N4  Palik Optical constants I';
if isempty(CCC_Si3N4) ; load retindice CCC_Si3N4;end;AAA=CCC_Si3N4;
%------------------------------------------------------
case 28.31;% indice de MgF2 (tetragonal ?) - ordinary ray - from Palik
bornes=[0.0459 ,2000]; ref='MgF2_no Palik Optical Constants : MgF2 no ko';    
if isempty(CCC_MgF2_ordinary) ; load retindice CCC_MgF2_ordinary;end;AAA=CCC_MgF2_ordinary;
%------------------------------------------------------
case 28.32;% indice de MgF2 (tetragonal ?) - extraordinary ray - from Palik
bornes=[0.0459, 2000]; ref='MgF2_ne Palik Optical Constants : MgF2 ne ke';  
if isempty(CCC_MgF2_extraordinary) ; load retindice CCC_MgF2_extraordinary;end;AAA=CCC_MgF2_extraordinary;
%------------------------------------------------------
case 4.13; % indice de GaSb (cubic) - from Palik, data added between 0.6 et 1 ï¿½m with available data from SOPRA Database (perticularly for k values)
bornes=[0.2067, 100]; ref='GaSb Palik Optical Constants'; 
if isempty(CCC_GaSb) ; load retindice CCC_GaSb;end;AAA=CCC_GaSb;
%------------------------------------------------------
case 28.41; % TiO2 ordinary - from Palik
bornes=[0.103, 13.2]; ref='TiO2_ordinary Palik optical constants : TiO2 ordinary (antiparallel)';
if isempty(CCC_TiO2_ordinary) ; load retindice CCC_TiO2_ordinary;end;AAA=CCC_TiO2_ordinary;
%------------------------------------------------------
case 28.42; % TiO2 extraordinary - from Palik
bornes=[0.103 ,13.2]; ref='TiO2_extraordinary Palik optical constants : TiO2 extraordinary (parallel)';
if isempty(CCC_TiO2_extraordinary) ; load retindice CCC_TiO2_extraordinary;end;AAA=CCC_TiO2_extraordinary;
%------------------------------------------------------
case 28.43; % TiO2 SOPRA 1 
bornes=[0.18, 1.5]; ref='TiO2  SOPRA  1';
if isempty(CCC_TiO2_SOPRA_1) ; load retindice CCC_TiO2_SOPRA_1;end;AAA=CCC_TiO2_SOPRA_1;
%------------------------------------------------------
case 28.44;  % TiO2 SOPRA 2
bornes=[0.22, 2.42]; ref='TiO2  SOPRA  2';
if isempty(CCC_TiO2_SOPRA_2) ; load retindice CCC_TiO2_SOPRA_2;end;AAA=CCC_TiO2_SOPRA_2;
%------------------------------------------------------
case 4.91;   % GaAs dope Si
bornes=[0.3, 3.0036]; ref='GaAs_dope_Si (n-type) J. Electromagnetic Analysis & Applications, 2010, 2, 357-361';
if isempty(CCC_GaAs_dope_Si) ; load retindice CCC_GaAs_dope_Si;end;AAA=CCC_GaAs_dope_Si;
%------------------------------------------------------
case 4.701; %  GaAs
bornes=[0.007999, 1000]; ref='GaAs Palik Optical Constants';
if isempty(CCC_GaAs) ; load retindice CCC_GaAs;end;AAA=CCC_GaAs;
%------------------------------------------------------
case 4.601;% Al20Ga80As
bornes=[0.496329588014981, 10]; ref='Al20Ga80As : S. Adachi : Properties of Gallium Arsenide {Emis Datareviews Series No 7, London : INSPEC, 1993}';
if isempty(CCC_Al20Ga80As) ; load retindice CCC_Al20Ga80As;end;AAA=CCC_Al20Ga80As;
%------------------------------------------------------
case 4.602; % Al315Ga685As
bornes=[0.496329588014981,2.48164794007491]; ref='Al315Ga685As : S. Adachi : Properties of Gallium Arsenide {Emis Datareviews Series No 7, London : INSPEC, 1993}';
if isempty(CCC_Al315Ga685As) ; load retindice CCC_Al315Ga685As;end;AAA=CCC_Al315Ga685As;
%------------------------------------------------------
case 4.6021; % Al30Ga70As
bornes=[0.000124,100]; ref='Al30Ga70As  Palik optical constants : Al30Ga70As recollé';
if isempty(CCC_Al30Ga70As) ; load retindice CCC_Al30Ga70As;end;AAA=CCC_Al30Ga70As;
%------------------------------------------------------
case 4.603;  % Al50Ga50As
bornes=[0.496329588014981,2.48164794007491]; ref='Al50Ga50As : S. Adachi : Properties of Gallium Arsenide {Emis Datareviews Series No 7, London : INSPEC, 1993}';
if isempty(CCC_Al50Ga50As) ; load retindice CCC_Al50Ga50As;end;AAA=CCC_Al50Ga50As;
%------------------------------------------------------
case 4.604;  % Al70Ga30As
bornes=[0.496329588014981,2.48164794007491]; ref='Al70Ga30As : S. Adachi : Properties of Gallium Arsenide {Emis Datareviews Series No 7, London : INSPEC, 1993}';
if isempty(CCC_Al70Ga30As) ; load retindice CCC_Al70Ga30As;end;AAA=CCC_Al70Ga30As;
%------------------------------------------------------
case 4.605;  %  Al80Ga20As
bornes=[0.2066,2.48]; ref='Al80Ga20As Palik Optical Constants (deux tables mises bout à bout) ';
if isempty(CCC_Al80Ga20As) ; load retindice CCC_Al80Ga20As;end;AAA=CCC_Al80Ga20As;
%------------------------------------------------------
case 4.606;  % AlAs
bornes=[0.496329588014981,2.48164794007491]; ref='AlAs : S. Adachi : Properties of Gallium Arsenide {Emis Datareviews Series No 7, London : INSPEC, 1993}';
if isempty(CCC_AlAs) ; load retindice CCC_AlAs;end;AAA=CCC_AlAs;

    
otherwise
error(['Milieu "' milieu '" inconnu!!'])
end  % switch parm
%------------------------------------------------------


lld=AAA(:,1);nn=AAA(:,2)+i*AAA(:,3);
if isempty(lambda);indice=[];else;indice=interp1(lld,nn,lambda,method,NaN*(1+1i));end; 
ldr=lld;ldi=lld;nr=real(nn);ni=imag(nn);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [n,bornes,ref,ldr,nr,ldi,ni]=indice_Collin_interpole(ld,parm,method);

persistent CCC_ld_nAg_nAu_nCr_nAsGa
n=ones(size(ld));lld=[];nn=[];

%------------------------------------------------------
% Donnees interpolees par S Collin
if isempty(CCC_ld_nAg_nAu_nCr_nAsGa) ; load retindice CCC_ld_nAg_nAu_nCr_nAsGa;end;ld_nAg_nAu_nCr_nAsGa=CCC_ld_nAg_nAu_nCr_nAsGa;
lld=ld_nAg_nAu_nCr_nAsGa(:,1);


switch parm;
case 101.7;nn=ld_nAg_nAu_nCr_nAsGa(:,2)+1i*ld_nAg_nAu_nCr_nAsGa(:,3);ref='Ag';
case 102.7;nn=ld_nAg_nAu_nCr_nAsGa(:,4)+1i*ld_nAg_nAu_nCr_nAsGa(:,5);ref='Au';
case 103.7;nn=ld_nAg_nAu_nCr_nAsGa(:,6)+1i*ld_nAg_nAu_nCr_nAsGa(:,7);ref='Cr';
case 104.7;nn=ld_nAg_nAu_nCr_nAsGa(:,8)+1i*ld_nAg_nAu_nCr_nAsGa(:,9);ref='GaAs';
case 114;nn=ld_nAg_nAu_nCr_nAsGa(:,10)+1i*ld_nAg_nAu_nCr_nAsGa(:,11);ref='ZnO';
case 114.1;nn=ld_nAg_nAu_nCr_nAsGa(:,12)+1i*ld_nAg_nAu_nCr_nAsGa(:,13);ref='ZnO_Al';
case 115;nn=ld_nAg_nAu_nCr_nAsGa(:,14)+1i*ld_nAg_nAu_nCr_nAsGa(:,15);ref='CdS';
case 116;nn=ld_nAg_nAu_nCr_nAsGa(:,16)+1i*ld_nAg_nAu_nCr_nAsGa(:,17);ref='Mo';
case 117;nn=ld_nAg_nAu_nCr_nAsGa(:,18)+1i*ld_nAg_nAu_nCr_nAsGa(:,19);ref='CIGS23';
case 117.1;nn=ld_nAg_nAu_nCr_nAsGa(:,20)+1i*ld_nAg_nAu_nCr_nAsGa(:,21);ref='CIGS31';
case 118;nn=ld_nAg_nAu_nCr_nAsGa(:,22)+1i*ld_nAg_nAu_nCr_nAsGa(:,23);ref='microSi';
case 118.1;nn=ld_nAg_nAu_nCr_nAsGa(:,24)+1i*ld_nAg_nAu_nCr_nAsGa(:,25);ref='a_Si';
case 119;nn=ld_nAg_nAu_nCr_nAsGa(:,26)+1i*ld_nAg_nAu_nCr_nAsGa(:,27);ref='CustomAl';
case 120;nn=ld_nAg_nAu_nCr_nAsGa(:,28)+1i*ld_nAg_nAu_nCr_nAsGa(:,29);ref='InP';
case 114.3;nn=ld_nAg_nAu_nCr_nAsGa(:,30)+1i*ld_nAg_nAu_nCr_nAsGa(:,31);ref='ZnO:Al2'; 
case 114.2;nn=ld_nAg_nAu_nCr_nAsGa(:,32)+1i*ld_nAg_nAu_nCr_nAsGa(:,33);ref='ZnOi2';
case 104.8;nn=ld_nAg_nAu_nCr_nAsGa(:,34)+1i*ld_nAg_nAu_nCr_nAsGa(:,35);ref='GaAs_C';
case 115.1;nn=ld_nAg_nAu_nCr_nAsGa(:,36)+1i*ld_nAg_nAu_nCr_nAsGa(:,37);ref='CdS2'; 
case 121;nn=ld_nAg_nAu_nCr_nAsGa(:,38)+1i*ld_nAg_nAu_nCr_nAsGa(:,39);ref='ITO';
case 118.2;nn=ld_nAg_nAu_nCr_nAsGa(:,40)+1i*ld_nAg_nAu_nCr_nAsGa(:,41);ref='a_Si_i';
case 118.3;nn=ld_nAg_nAu_nCr_nAsGa(:,42)+1i*ld_nAg_nAu_nCr_nAsGa(:,43);ref='a_SiC_p';
case 122;nn=ld_nAg_nAu_nCr_nAsGa(:,44)+1i*ld_nAg_nAu_nCr_nAsGa(:,45);ref='ITO_2';   
   
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f=find(~isnan(nn));lld=lld(f);nn=nn(f);bornes=[min(lld),max(lld)];
if isempty(ld);n=[];else; n=interp1(lld,nn,ld,method,NaN*(1+1i));end;
ref=[ref,'  Donnees interpollees par S Collin corrigées JPH'];
ldr=lld;ldi=lld;nr=real(nn);ni=imag(nn);
    





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [indice,bornes,ref,ldr,nr,ldi,ni]=indice_mathieu(lambda,parm,method);
% Transmis par Mathieu 2014

% {13.24,13.241,26,26.1,4.15,4.81,3.26,22.1,27,23.1,23.11,30,31,32,33,40,41,42,43,44}
persistent MMM_TUNGSTEN_24 MMM_TUNGSTEN_241 MMM_MOLYBDENE	MMM_MOLYBDENE_1	MMM_Si3N4	MMM_TiO2 ...
MMM_CHROME	MMM_MgO	MMM_Tantalum	MMM_SiC_Mesure_ellipso	MMM_SiC_Mesure_ellipso_fite	MMM_Glass_HEF	MMM_Tungsten_HEF	MMM_MgO_HEF ...
MMM_TiO2_HEF MMM_MgO_Material_Design	MMM_SiC_cubic_Material_Design	MMM_TiO2_rutile_Material_Design	MMM_W_Material_Design	MMM_WO2_Material_Design;
indice=ones(size(lambda));lld=[];nn=[];


switch parm;
% -----------------------------------------
case 13.24;
bornes=[.01     24.80];ref='W   Palick  et Electronic Handbook of Optical Constant of solid version 1.0 HOC';
if isempty(MMM_TUNGSTEN_24) ; load retindice MMM_TUNGSTEN_24;end;AAA=MMM_TUNGSTEN_24;
% -----------------------------------------
case 13.241;
bornes=[.01     4.5];ref='W   Palick  et Electronic Handbook of Optical Constant of solid version 1.0 HOC';
if isempty(MMM_TUNGSTEN_241) ; load retindice MMM_TUNGSTEN_241;end;AAA=MMM_TUNGSTEN_241;
% -----------------------------------------
case 26;
bornes=[0.001   12.40];ref='Mo   Palick  et Electronic Handbook of Optical Constant of solid version 1.0 HOC';
if isempty(MMM_MOLYBDENE) ; load retindice MMM_MOLYBDENE;end;AAA=MMM_MOLYBDENE;
% -----------------------------------------
case 26.1;
bornes=[0.001   4.5];ref='Mo   Palick  et Electronic Handbook of Optical Constant of solid version 1.0 HOC';
if isempty(MMM_MOLYBDENE_1) ; load retindice MMM_MOLYBDENE_1;end;AAA=MMM_MOLYBDENE_1;
% -----------------------------------------
case 4.15;
bornes=[0.200   6];ref='Si3N4 (Palik ? Mathieu)';
if isempty(MMM_Si3N4) ; load retindice MMM_Si3N4;end;AAA=MMM_Si3N4;
% -----------------------------------------
case 4.81;
bornes=[0.001   4.5];ref='TiO2   Palick  et Electronic Handbook of Optical Constant of solid version 1.0 HOC';
if isempty(MMM_TiO2) ; load retindice MMM_TiO2;end;AAA=MMM_TiO2;
% -----------------------------------------
case 3.26;
bornes=[0.01   20];ref='Cr   Palick  et Electronic Handbook of Optical Constant of solid version 1.0 HOC';
if isempty(MMM_CHROME) ; load retindice MMM_CHROME;end;AAA=MMM_CHROME;
% -----------------------------------------
case 22.1;
bornes=[0.01653   20];ref='MgO   Palick  et Electronic Handbook of Optical Constant of solid version 1.0 HOC';
if isempty(MMM_MgO) ; load retindice MMM_MgO;end;AAA=MMM_MgO;
% -----------------------------------------
case 27;
bornes=[0.000124   125];ref='Ta   Palick  et Electronic Handbook of Optical Constant of solid version 1.0 HOC';
if isempty(MMM_Tantalum) ; load retindice MMM_Tantalum;end;AAA=MMM_Tantalum;
% -----------------------------------------
case 23.12;
bornes=[0.262   10];ref='SiC Mesure ellipso';
if isempty(MMM_SiC_Mesure_ellipso) ; load retindice MMM_SiC_Mesure_ellipso;end;AAA=MMM_SiC_Mesure_ellipso;
% -----------------------------------------
case 23.11;
bornes=[0.262   10];ref='SiC Mesure ellipso fitée';
if isempty(MMM_SiC_Mesure_ellipso_fite) ; load retindice MMM_SiC_Mesure_ellipso_fite;end;AAA=MMM_SiC_Mesure_ellipso_fite;
% -----------------------------------------
case 30;
bornes=[0.300   2];ref='Verre HEF, mesures R&T, fit avec SCOUT';
if isempty(MMM_Glass_HEF) ; load retindice MMM_Glass_HEF;end;AAA=MMM_Glass_HEF;
% -----------------------------------------
case 31;
bornes=[0.300   2];ref='W HEF, mesures R&T, fit avec SCOUT';
if isempty(MMM_Tungsten_HEF) ; load retindice MMM_Tungsten_HEF;end;AAA=MMM_Tungsten_HEF;
% -----------------------------------------
case 32;
bornes=[0.300   2];ref='MgO HEF, mesures R&T, fit avec SCOUT';
if isempty(MMM_MgO_HEF) ; load retindice MMM_MgO_HEF;end;AAA=MMM_MgO_HEF;
% -----------------------------------------
case 33;
bornes=[0.300   2];ref='TiO2 HEF, mesures R&T, fit avec SCOUT';
if isempty(MMM_TiO2_HEF) ; load retindice MMM_TiO2_HEF;end;AAA=MMM_TiO2_HEF;
% -----------------------------------------
case 40;
bornes=[0.100   22];ref='MgO Material Design';
if isempty(MMM_MgO_Material_Design) ; load retindice MMM_MgO_Material_Design;end;AAA=MMM_MgO_Material_Design;
% -----------------------------------------
case 41;
bornes=[0.200   13];ref='SiC cubic Material Design';
if isempty(MMM_SiC_cubic_Material_Design) ; load retindice MMM_SiC_cubic_Material_Design;end;AAA=MMM_SiC_cubic_Material_Design;
% -----------------------------------------
case 42;
bornes=[0.200   13];ref='TiO2 rutile Material Design';
if isempty(MMM_TiO2_rutile_Material_Design) ; load retindice MMM_TiO2_rutile_Material_Design;end;AAA=MMM_TiO2_rutile_Material_Design;
% -----------------------------------------
case 43;
bornes=[0.200   6];ref='W Material Design';
if isempty(MMM_W_Material_Design) ; load retindice MMM_W_Material_Design;end;AAA=MMM_W_Material_Design;
% -----------------------------------------
case 44;
bornes=[0.200   6];ref='WO2 Material Design';
if isempty(MMM_WO2_Material_Design) ; load retindice MMM_WO2_Material_Design;end;AAA=MMM_WO2_Material_Design;

otherwise
error(['Milieu  inconnu!!'])
end  % switch parm
%------------------------------------------------------


lld=AAA(:,1);nn=AAA(:,2)+i*AAA(:,3);
if isempty(lambda);indice=[];else;indice=interp1(lld,nn,lambda,method,NaN*(1+1i));end; 
ldr=lld;ldi=lld;nr=real(nn);ni=imag(nn);






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [indice,bornes,ref,ldr,nr,ldi,ni]=indiceor(ld);
lld=[];nn=[];
bornes=[.55,.95];ref='Au  approximation polynomiale de degre 4 de Corn   (2.04) ';
% ld microns 
indice=(29.0387-140.431*ld+254.394*ld.^2-203.936*ld.^3+61.1781*ld.^4)-i*(5*ld.^2-16.3*ld+4.9253);
indice(ld<bornes(1)|ld>bornes(2))=nan+i*nan;
ldr=lld;ldi=lld;nr=real(nn);ni=imag(nn);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [indice,bornes,ref,ldr,nr,ldi,ni]=indiceor_drude(ld);
lld=[];nn=[];
bornes=[-inf,inf];ref='Au TRES MAUVAIS *** modele de Drude  benchmark reseau  (mauvaise approximation)';
% ld microns
c=2.99792458e14;omegap=1.374e16/(2*pi*c);gama=3.21e14/(2*pi*c);
indice=retsqrt(1-omegap^2*(ld.^2./(1+i*gama*ld)),-1);
ldr=lld;ldi=lld;nr=real(nn);ni=imag(nn);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [indice,bornes,ref,ldr,nr,ldi,ni]=indiceor_drude_lorentz(ld);
lld=[];nn=[];
bornes=[-inf,inf];ref='Au modele de Drude Lorentz: F.Kaminski,V.Sandoghdar,and M.Agio Journal of Computional and Theoretical Nanoscience Vol 4 635-643 2007';
% ld microns
indice=retsqrt(8.0969-ld.^2*49.0368./(1+0.0219i*ld)-3.4233*ld.^2./(1-ld.^2*4.6669+0.5021i*ld),-1);
%indice=retsqrt(8.0969-ld.^2*49./(1+0.022i*ld)-3.4233*ld.^2./(1-ld.^2*4.667+0.5i*ld),-1);
ldr=lld;ldi=lld;nr=real(nn);ni=imag(nn);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [indice,bornes,ref,ldr,nr,ldi,ni]=indiceSiC_drude(ld);
lld=[];nn=[];
bornes=[-inf,inf];ref='SiC  modele de Lorentz donnees Palik transmis par:J le Gall M Olivier JJ Greffet Phys Rev B vol 55 number 15 15 april 1997'; 
omegaL=969;omegaT=793;Gama=-4.76;eps_inf=6.7;
indice=retsqrt(eps_inf*(1+(omegaL^2-omegaT^2)./(omegaT^2-(10000./ld).^2+i*Gama*(10000./ld))),-1);
ldr=lld;ldi=lld;nr=real(nn);ni=imag(nn);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [indice,bornes,ref,ldr,nr,ldi,ni] = indiceSi_Lorentz(ld,parm);% transmis par Mathieu remis en -i omega t par JPH
lld=[];nn=[];
bornes=[-inf,inf];
% ld microns
c=3e8;
switch parm;
case  23.2; 
ref='SiC_np Double oscillateur de Lorentz';
f1 = 4.0179362; wo1 = 0.0989483; g1 = 0.0029082;
f2 = 1.0499011; wo2 = 0.7416009; g2 = 0.4762199;
case 23.3;
ref='SiC_non_poli Double oscillateur de Lorentz';
f1 = 3.44721329 ; wo1 = 0.0985046; g1 = 0.0007774;
f2 = 0.331191; wo2 = 0.2799479; g2 = 0.2521959;
end;
eps_inf = 6.84;
w = c./(ld*1e-6);
w = w.*4.1357e-15;
eps = eps_inf +(f1*wo1^2)./(wo1^2-w.^2-1i*g1.*w) + (f2*wo2^2)./(wo2^2-w.^2-1i*g2.*w);
indice=retsqrt(eps,-1);
ldr=lld;ldi=lld;nr=real(nn);ni=imag(nn);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [indice,bornes,ref,ldr,nr,ldi,ni]=indiceSiO2_Lorentz(ld);
lld=[];nn=[];
bornes=[-inf,inf];ref='SiO2 modele de Lorentz d''apres 4.41 fité par JPH valable dans le domaine  1.5385, 20';
Omega_L=1220.8;
Omega_T=1048.7;
Gamma=71.3901;
eps_inf=2.0948;
indice=retsqrt(eps_inf*(1+(Omega_L^2-Omega_T^2)./((Omega_T-(1.e4./ld)).*(Omega_T+(1.e4./ld))-1i*Gamma*(1.e4./ld))));
ldr=lld;ldi=lld;nr=real(nn);ni=imag(nn);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [indice,bornes,ref,ldr,nr,ldi,ni]=indiceag_drude_Alex(ld);
lld=[];nn=[];
bornes=[-inf,inf];ref='Ag modele de Drude donné par A Archambault,F Marquier, JJ Greffet PRB 82 ,035411 (2010)'; 
% ld microns
eV=retdb(ld,'ld_2_eV');
indice=retsqrt(5-9.1^2./(eV.^2+.021i*eV)+1.8i./eV,-1);
ldr=lld;ldi=lld;nr=real(nn);ni=imag(nn);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [indice,bornes,ref,ldr,nr,ldi,ni]=indiceag_drude_Johnson(ld);
lld=[];nn=[];
bornes=[-inf,inf];ref='Ag  modéle de Drude donné par P.B. Johnson. R.W. Christy transmis par Christophe Sauvan '; 
% ld microns
indice=retsqrt(1-ld.^2./(.138^2*(1+.017i*ld)),-1);
ldr=lld;ldi=lld;nr=real(nn);ni=imag(nn);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [indice,bornes,ref,ldr,nr,ldi,ni]=indiceag_drude(ld);
lld=[];nn=[];
bornes=[.75,.9];ref='Ag  ASSEZ MEDIOCRE *** modele de Drude ''comment'' G Gay,Alloschery,Weiner .Nature Phys (assez mediocre ...) '; 
% ld microns
c=2.99792458e14;omegap=1.3552e16/(2*pi*c);gama=1.9944e14/(2*pi*c);
indice=retsqrt(3.2938-omegap^2*(ld.^2./(1+i*gama*ld)),-1);
indice(ld<bornes(1)|ld>bornes(2))=nan+i*nan;
ldr=lld;ldi=lld;nr=real(nn);ni=imag(nn);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [indice,bornes,ref,ldr,nr,ldi,ni]=indiceag_drude_benchmark(ld);
lld=[];nn=[];
bornes=[.4,2];ref='Ag  TRES MAUVAIS *** modele de Drude ''benchmark''  ';
% ld microns
c=3e14;omegap=1.3388e16/(2*pi*c);gama=.707592e14/(2*pi*c);
indice=retsqrt(3.36174-omegap^2*(ld.^2./(1+i*gama*ld)),-1);
indice(real(ld)<bornes(1)|real(ld)>bornes(2))=nan+i*nan;
ldr=lld;ldi=lld;nr=real(nn);ni=imag(nn);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [indice,bornes,ref,ldr,nr,ldi,ni]=indicealu_drude(ld);
lld=[];nn=[];
bornes=[-inf,inf];ref='Al modele de Drude PhysRevB.68.245405';
% ld microns
c=299792458;omegap=1.747e16/(c*1.e6);gama=7.596e13/(c*1.e6);% multiplier par c*1.e6=2.997924580000000e+014 pour avoir omegap,gama en rd/s
indice=retsqrt(1-omegap.^2./((2*pi./ld).*((2*pi./ld)+1i*gama)),-1);
ldr=lld;ldi=lld;nr=real(nn);ni=imag(nn);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [indice,bornes,ref,ldr,nr,ldi,ni]=indice_Si(ld,T,cas);  
if ischar(T);T=293;end;
switch cas;
case 0; %silicium    
bornes=[1.1,5.6];
S1=10.4907-2.0802e-4*T+4.21694e-6*T^2-5.82298e-9*T^3+3.44688e-12*T^4;
S2=-1346.61+29.1664*T-.278724*T^2+1.05939e-3*T^3-1.35089e-6*T^4;
S3=4.42827e7-1.76213e6*T-7.61575e4*T^2+678.414*T^3+103.243*T^4;
ld1=.299713-1.14234e-5*T+1.67134e-7*T^2-2.51049e-10*T^3+2.32484e-14*T^4;
ld2=-3.5171e3+42.3892*T-.357957*T^2+1.17504e-3*T^3-1.13212e-6*T^4;
ld3=1.714e6-1.44984e5*T-6.90744e3*T^2-39.3699*T^3+23.577*T^4;
ref='Si Temperature dependent refractive index of silicon and germanium Frey Leviton Madison,Proc,SPIE Vol,6273,6272J'; 

case 1; %germanium    
bornes=[1.9,5.5];
S1=13.9723+2.52809e-3*T-5.02195e-6*T^2+2.22604e-8*T^3-4.86238e-12*T^4;
S2=.452096-3.09197e-3*T+2.16895e-5*T^2-6.02290e-8*T^3+4.12038e-11*T^4;
S3=751.447-14.2843*T-.238093*T^2+2.96047e-3*T^3-7.73454e-6*T^4;
ld1=.386367+2.01871e-4*T-5.93448e-7*T^2-2.27923e-10*T^3+5.37423e-12*T^4;
ld2=1.08843+1.16510e-3*T-4.97284e-6*T^2+1.12357e-8*T^3+9.40201e-12*T^4;
ld3=-2893.19-.967948*T-.527016*T^2+6.49364e-3*T^3-1.95162e-5*T^4;
end;
indice=retsqrt(1+S1*ld.^2./(ld.^2-ld1^2)+S2*ld.^2./(ld.^2-ld2^2)+S3*ld.^2./(ld.^2-ld3^2),-1);

ref='Ge Temperature dependent refractive index of silicon and germanium Frey Leviton Madison,Proc,SPIE Vol,6273,6272J'; 
indice(ld<bornes(1)|ld>bornes(2))=nan+i*nan;
[lld,nn,ldr,ldi,nr,ni]=deal([]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [indice,bornes,ref,ldr,nr,ldi,ni]=indice_algaas(cas,lambda,x)
% Cette fonction retourne l'indice du matériau AlxGa1-xAs en fonction du paramètre x
% et de la longueur d'onde exprimée en microns
% Référence : "Optoélectronique", Emmanuel Rosencher et Borge Vinter, Masson, chapitre 12.7, P332
ref='Optoélectronique, Emmanuel Rosencher et Borge Vinter, Masson, chapitre 12.7, P332';
bornes=[-inf,inf];
E = 1.24./lambda; % E en eV, lambda en microns
[lld,nn,ldr,ldi,nr,ni]=deal([]);
if ischar(x);x=0;end;

switch cas;
case 0;ref=['AlxGA1-xAs ',ref];
Eoe = 3.65 + 0.871*x + 0.179*(x^2);
Ed  = 36.1 - 2.45*x;
Eg  = 1.424 + 1.266*x + 0.26*(x^2);
case 1;ref=['PxAs1-xGa ',ref];
Eoe = 3.65 + 0.721*x + 0.139*(x^2);
Ed  = 36.1+0.35*x;
Eg = 1.441 + 1.091*x + 0.21*(x^2);
case 2;ref=['GaxIn1-xP ',ref];
Eoe = 3.391 + 0.524*x + 0.595*(x^2);
Ed = 28.91 + 7.54*x;
Eg = 1.34 + 0.668*x + 0.758*(x^2);
end;
Ef = sqrt((2*(Eoe^2)) - Eg^2);
eta = (pi/2) * (Ed/((Eoe^3)*(Eoe^2 - Eg^2)));
M_1 = (eta/(2*pi))*(Ef^4 - Eg^4);
M_3 = (eta/pi) * (Ef^2 - Eg^2);
epsilon = 1 + M_1 + M_3*(E.^2) + ((eta/pi)*log((Ef^2-E.^2)./(Eg^2 - E.^2))).*E.^4; % Formule d'Afromowitz
indice = retsqrt(epsilon,-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [indice,bornes,ref,ldr,nr,ldi,ni]=indice_GaAs_Simon(ld,Nd)
% Référence : "Optoélectronique", Emmanuel Rosencher et Borge Vinter, Masson, chapitre 12.7, P332
% Nd en densité electronique en nb par microns cube
% 1.5e18 /cm^3 --> 1.5e6 /mu^3
% retindice(35,4.63) GaAs pas dopé
% retindice(35,4.63,1.5e6) GaAs dopé Nd=1.5e6 /mu^3=1.5e18 /cm^3 Gamma_e=60
% retindice(35,4.63,[1.5e6,40]) GaAs dopé Nd=1.5e6 /mu^3=1.5e18 /cm^3 Gamma_e=40

ref='GaAs  formule de  Simon Vassant Palik  Electronic Handbook of Optical Constant of solid ISBN 0-12-544420-60 p429';
bornes=[-inf,inf];
[lld,nn,ldr,ldi,nr,ni]=deal([]);
if ischar(Nd);Nd=0;end;
% definition des indices (parametre de Palik)
Omega_L=292.1;
Omega_T=268.7;
Gamma=2.4;
eps_inf=11;
if length(Nd)==1;Gamma_e=60;else;Gamma_e=Nd(2);Nd=Nd(1);end;
Nd=Nd*1.e18;%passage en m-3  %1mu3 = 1e-18m3
%Omega_p=sqrt(Nd*e^2/(eps0*m_star))*1.e-2/(2*pi*c);% pour passer en cm-1
%c=2.99792458e8;e=1.60219e-19;eps0=8.85419e-12;m_star=0.067*9.10953e-31;Omega_p2_s_Nd=(e^2/(eps0*m_star))   *1.e-4/(2*pi*c).^2;% pour passer en cm-1
Omega_p2_s_Nd=1.338774034689134e-018;

Omega_Z=roots([-1,-1i*Gamma,Omega_T^2]);

Omega=1.e4./ld;
%indice=retsqrt(eps_inf*(1+(Omega_L^2-Omega_T^2)./((Omega_T-Omega).*(Omega_T+Omega)-1i*Gamma*Omega))-Omega_p2./(Omega.^2+1i*Gamma_e*Omega),-1);
indice=retsqrt(eps_inf*(1-(Omega_L^2-Omega_T^2)./( (Omega-Omega_Z(1)).*(Omega-Omega_Z(2)) ))-(Nd*Omega_p2_s_Nd)./(Omega.*(Omega+1i*Gamma_e)),-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [indice,bornes,ref,ldr,nr,ldi,ni]=indice_drudePBA(ld,parm);
lld=[];nn=[];
bornes=[-inf,inf];ref='modele de Drude Ordal transmis par P Ben-Abdallah Optical properties of fourteen metals APPLIED OPTIC 1985 pp4493-4499';
switch parm;
case 1.6;Omegap=13.69;gamma=2.73;ref=['Ag TRES MAUVAIS ***',ref];
case 1.61;Omegap=11.381;gamma=9.6338;ref=['Ag ',ref,' correction JPH 2 2014'];
case 10.6;Omegap=22.42;gamma=12.43;ref=['Al TRES MAUVAIS ***',ref];
case 2.6;Omegap=13.71;gamma=4.05;ref=['Au TRES MAUVAIS ***',ref];
case 2.61;Omegap=10.965;gamma=8.5696;ref=['Au ',ref,' correction JPH 2 2014'];
case 15.6;Omegap=6.03;gamma=5.56;ref=['Co  TRES MAUVAIS ***',ref];
case 12.6;Omegap=11.23;gamma=1.38;ref=['Cu TRES MAUVAIS ***',ref];
case 12.61;Omegap=10.359;gamma=9.6937;ref=['Cu ',ref,' correction JPH 2 2014'];
case {11.6,18.6};Omegap=3.82;gamma=7.20;ref=['Ti  TRES MAUVAIS ***',ref];
case 18.61;Omegap=3.9294;gamma=4.4356;ref=['Ti  ',ref,' correction JPH 2 2014'];
end;
c=retconstantes('c');
%Omegap=1.e9*Omegap/c;gamma=1.e7*gamma/c;
Omegap=1.e9*Omegap/c;gamma=1.e7*gamma/c;
indice=retsqrt(1-Omegap^2./((2*pi./ld).*((2*pi./ld)+1i*gamma)),-1);
ldr=lld;ldi=lld;nr=real(nn);ni=imag(nn);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [indice,bornes,ref,ldr,nr,ldi,ni]=indice_Rakic_metaux(ld,cas);
lld=[];nn=[];
bornes=[-inf,inf];ref='A.D.Rakic,A B.Djurisic,J.M.Elazar,M.L.Majewski,''Optical properties of metallic films for vertical-cavity optoelectronic devices,'' Appl. Opt. 37, 5271-5283 (1998)';
ldr=lld;ldi=lld;nr=real(nn);ni=imag(nn);
% cette fonction a ete teste avec le programme en Pyton ldbb.
% A Python script written by Aaron Webster for calculating optical constants of several metals using the
% Drude, Lorentz-Drude, or Brendel-Bormann models.
%
% 1:Ag  2:Au  3:Cu  4:Al  5:Be  6:Cr  7:Ni  8:Pd  9:Pt  10:Ti  11: W
eV=retdb(ld,'ld_2_eV');
switch(floor(cas));
case 1;num=1;ref=['Ag ',ref];	
case 2;num=2;ref=['Au ',ref];	
case 12;num=3;ref=['Cu ',ref];	
case 10;num=4;ref=['Al ',ref];	
case 25;num=5;ref=['Be ',ref];	
case 3;num=6;ref=['Cr ',ref];	
case 16;num=7;ref=['Ni ',ref];	
case 17;num=8;ref=['Pd ',ref];	
case 24;num=9;ref=['Pt ',ref];	
case 18;num=10;ref=['Ti ',ref];	
case 13;num=11;ref=['W ',ref];	
end;	
Omegap={9.01,9.03,10.83,14.98,18.51,10.75,15.92,9.72,9.59,7.29,13.22};
	
Modele=round(100*(cas-floor(cas)-.9));	
if Modele==2; % BB
f0={.821,.770,.562,.526,.081,.154,.083,.33,.333,.126,.197};
Gama0={.049,.050,.03,.047,.035,.048,.022,.009,.08,.067,.057};
f={[.05,.133,.051,.467,4],...
	[.054,.05,.312,.719,1.648],...
	[.076,.081,.324,.726],...
	[.213,.06,.182,.014],...
	[.066,.067,.346,.311],...
	[.338,.261,.817,.105],...
	[.357,.039,.127,.654],...
	[.769,.093,.309,.409],...
	[.186,.665,.551,2.214],...
	[.427,.218,.513,.0002],...
	[.006,.022,.136,2.648]};
Gama={[.189,.067,.019,.117,.052],...
	[.074,.035,.083,.125,.179],...
	[.056,.047,.113,.172],...
	[.312,.315,1.587,2.145,],...
	[2.956,3.962,2.398,3.904],...
	[4.256,3.957,2.218,6.983],...
	[2.82,.12,1.822,6.637,],...
	[2.343,.497,2.022,.119,],...
	[.498,1.851,2.604,2.891],...
	[1.877,.1,.615,4.109],...
	[3.689,.277,1.433,4.555]};

Omega={[2.025,5.185,4.343,9.809,18.56],...
	[.218,2.885,4.069,6.137,27.97],...
	[.416,2.849,4.819,8.136],...
	[.163,1.561,1.827,4.495],...
	[.131,.469,2.827,4.318],...
	[.281,.584,1.919,6.997],...
	[.317,1.059,4.583,8.825],...
	[.066,.502,2.432,5.987],...
	[.782,1.317,3.189,8.236],...
	[1.459,2.661,.805,19.86],...
	[.481,.985,1.962,5.442]};

Sigma={[1.894,.665,.189,1.170,.516],...
	[.742,.349,.830,1.246,1.795],...
	[.562,.469,1.131,1.719],...
	[.013,.042,.256,1.735],...
	[.277,3.167,1.446,.893],...
	[.115,.252,.225,4.903],...
	[.606,1.454,.379,.510],...
	[.694,.027,1.167,1.331],...
	[.031,.096,.766,1.146],...
	[.463,.506,.799,2.854],...
	[3.754,.059,.273,1.912]};

Omegap=Omegap{num};f0=f0{num};Gama0=Gama0{num};f=f{num};Gama=Gama{num};Omega=Omega{num};Sigma=Sigma{num};
ep=1-f0*Omegap^2./(eV.*(eV+1i*Gama0));
for ii=1:length(f);
A=retsqrt(eV.*(eV+1i*Gama(ii)),0);
A=retsqrt(eV.*(eV+1i*Gama(ii)),0);
ep=ep+(f(ii)*Omegap^2/(sqrt(8*pi)*Sigma(ii)))*(  retoc( (A-Omega(ii))/(sqrt(2)*Sigma(ii))) - retoc((-A-Omega(ii))/(sqrt(2)*Sigma(ii)))   )./A;	
% 
% 	Am=(A-Omega(ii))/(sqrt(2)*Sigma(ii));
% 	Ap=(A+Omega(ii))/(sqrt(2)*Sigma(ii));
% 	Xm=1i*Am*sqrt(2i/pi);[ym,yym,yyym,parmm]=retfresnel(Xm);
% 	Xp=1i*Ap*sqrt(2i/pi);[yp,yyp,yyyp,parmp]=retfresnel(Xp);
% 
% 	% Wm=exp(-Am.^2).*(1+sqrt(-2i)*yyym)-sqrt(-2i)*((i/pi)./Xm+.5i*yym./Xm);
% 	% Wp=exp(-Ap.^2).*(1+sqrt(-2i)*yyyp)-sqrt(-2i)*((i/pi)./Xp+.5i*yyp./Xp);
% 
% 	[fm,ffm]=retfind(parmm==1);
% 	[fp,ffp]=retfind(parmp==1);retsize(fm,ffm,fp,ffp)
% 	Wm=zeros(size(A));Wp=zeros(size(A));
% 	Wm(ffm)=exp(-Am(ffm).^2).*(1+sqrt(-2i)*yyym(ffm))-sqrt(-2i)*((i/pi)./Xm(ffm)+.5i*yym(ffm)./Xm(ffm));
% 	Wp(ffp)=exp(-Ap(ffp).^2).*(1+sqrt(-2i)*yyyp(ffp))-sqrt(-2i)*((i/pi)./Xp(ffp)+.5i*yyp(ffp)./Xp(ffp));
% 	Wm(fm)=exp(-Am(fm).^2).*(1+sqrt(-2i)*retfresnel(1i*Am(fm)*sqrt(2i/pi)));
% 	Wp(fp)=exp(-Ap(fp).^2).*(1+sqrt(-2i)*retfresnel(1i*Ap(fp)*sqrt(2i/pi)));
% 
% 	% reterfc=@(x) 1+sqrt(-2i)*retfresnel(-x*sqrt(2i/pi)) 	
% 	% Wm1=exp(-Am.^2).*(1+sqrt(-2i)*retfresnel(1i*Am*sqrt(2i/pi)));Wm1(isnan(Wm1))=0;
% 	% Wp1=exp(-Ap.^2).*(1+sqrt(-2i)*retfresnel(1i*Ap*sqrt(2i/pi)));Wp1(isnan(Wp1))=0;
% 	% erm=retcompare(Wm,Wm1),
% 	% erm=retcompare(Wp,Wp1),
% 	ep=ep+1i*sqrt(pi)*f(ii)*Omegap^2/(2*sqrt(2)*Sigma(ii))*(Wm+Wp)./A;
end;

else; % LD
f0={.845,.760,.575,.523,.084,.168,.096,.33,.333,.148,.206 };
Gama0={.048,.053,.030,.047,.035,.047,.048,.008,.08,.082,.064};
f={[.065,.124,.011,.840,5.646],...
	[.024,.01,.071,.601,4.384],...
	[.061,.104,.723,.638],...
	[.227, .05,.166,.03],...
	[.031,.140,.53,.130],...
	[.151,.150,1.149,.825],...
	[.1,.135,.106,.729],...
	[.649,.121,.638,.453],...
	[.191,.659,.547,3.576],...
	[.899,.393,.187,.001]...
	[.054,.166,.706,2.590]};
	
Gama={[3.886,.452,.065,.916,2.419],...
	[.241,.345,.870,2.494,2.214],...
	[.378,1.056,3.213,4.305],...
	[.333,.312,1.351,3.382],...
	[1.664,3.395,4.454,1.802],...
	[3.175,1.305,2.676,1.335],...
	[4.511,1.334,2.178,6.292],...
	[2.95,.555,4.621,3.236],...
	[.517,1.838,3.668,8.517],...
	[2.276,2.518,1.663,1.762],...
	[.53,1.281,3.332,5.836]};

Omega={[.816,4.481,8.185,9.083,20.29],...
	[.415,.83,2.969,4.304,13.32],...
	[.291,2.957,5.3,11.18],...
	[.162,1.544,1.808,3.473],...
	[.1,1.032,3.183,4.604],...
	[.121,.543,1.97,8.775],...
	[.174,.582,1.597,6.089],...
	[.336,.501,1.659,5.715],...
	[.78,1.314,3.141,9.249,],...
	[.777,1.545,2.509,19.43],...
	[1.004,1.917,3.58,7.498]};
Omegap=Omegap{num};f0=f0{num};Gama0=Gama0{num};f=f{num};Gama=Gama{num};Omega=Omega{num};
ep=1-f0*Omegap^2./(eV.*(eV+1i*Gama0));
for ii=1:length(f);
ep=ep+f(ii)*Omegap^2./(Omega(ii)^2-eV.^2+(-1i*Gama(ii))*eV);
end;
end;
indice=retsqrt(ep,-1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [indice,bornes,ref,ldr,nr,ldi,ni]=indice_Nordlander(ld,cas);
lld=[];nn=[];
bornes=[-inf,inf];ref='Modele de Lorentz à 4 termes Normander Effective dielectric function for FDTD... Chemical Physics Letters 446 (2007) 115-118 Transmis par JJG';
ldr=lld;ldi=lld;nr=real(nn);ni=imag(nn);
eV=retdb(ld,'ld_2_eV');

switch floor(cas);
case 2;ref=['Au ',ref];
Eps_inf=1;Sig=1355.01;A=[-8.577e4,-2.875,-997.6,-1.630];B=[-1.156e4,0,-3090,-4.409];C=[5.557e7,2.079e3,6.921e5,26.15];
case 1;ref=['Ag ',ref];
%Eps_inf=1;Sig=3157.56;A=[-1.160e5,-4.252,-0.496,-2.118];B=[-3050,-0.8385,-13.85,-10.23];C=[3.634e8,112.2,1.815,14.31];
Eps_inf=1;Sig=3158.56;3157.56;A=[-1.160e5,-4.252,-0.496,-2.118];B=[-3050,-0.8385,-13.85,-10.23];C=[3.634e8,112.2,1.815,14.31];
end;
ep=Eps_inf-Sig./(1i*eV);for ii=1:length(C);ep=ep+C(ii)./(eV.^2-A(ii)*1i*eV+B(ii));end;
indice=retsqrt(ep,-1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [n,bornes,ref]=indicefibre(ld,d,n2,n1);
bornes=[nan,nan];ref='Fibre formule approchee puis recherche par retmarcuse';
n2=real(n2);n1=real(n1);
% d et ld ont la meme longueur  n2 n1 reels n2/n1 peut varier entre 1.18 et 5
if isempty(ld)|isempty(d);n=[];return;end;
if length(ld)<length(d);ld=ld(1)*ones(size(d));end;
if length(d)<length(ld);d=d(1)*ones(size(ld));end;
if n2>5*n1;aa=[-7.3757  ;  5.4352 ;   2.3411  ;  1.1996];
% else;aa=[ -1.2058    5.4116  -11.8752;    0.7901   -3.2581    7.5089;   -0.1713    1.6521   -1.1545;   -0.1635    0.8319    0.2175];
else;aa=[-1.3744    6.4535  -12.8696;    0.8809   -3.8194    8.0447;   -0.1559    1.5568   -1.0635;   -0.2353    1.2757   -0.2061];
end;
a=zeros(1,size(aa,1));for ii=1:size(aa,1);a(ii)=polyval(aa(ii,:),n2/n1);end;
x=n2*d./ld;  % (n2/n1)*(d*n1/ld)
n=exp(-exp(a(1)*x+a(2))).* (1-1./(a(3).*x.^2+a(4)));
n=n1+(n2-n1)*n;
for ii=1:length(d);k0=2*pi/ld(ii);
if abs(n(ii)-n1)<10*eps;n(ii)=n1;
else;prv=real(retmarcuse([n1,n2],d(ii)/2,n(ii)*k0,1,[],[],k0)/k0);
if isnan(prv);if abs(n(ii)-n1)<abs(n(ii)-n2);n(ii)=n1;else;n(ii)=n2;end;else;n(ii)=prv;end;	
end;
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [n,bornes,ref]=indice_BQ(ld,ld0,ctho,n_ext,h_boite,r_boite,F);
bornes=[-inf,inf];ref='BQ formule BQ modele pour la polarisabilité de l électron élastiquement lié, + géométrie de la boite donnée par Florian ';
if nargin<7;F=1/3;end;
% h_boite=.004;r_boite=.01;ctho=3e14*1.3e-9;ld0=.91;F force de l'oscillateur (par defaut 3F=1)
%  A*Volume_BQ = charge_electron ^2  mu_0 / masse_electron = ((1.60219e-19)^2 *4.e-7*pi/9.10953e-31)* 1.e6 =3.541130460912364e-008     unite =longueur en microns
A=3*F*3.541130460912364e-008/(pi*r_boite^2*h_boite);
%n=retsqrt(n_ext^2+A./((2*pi./ld).^2-(2*pi./ld0).^2-1i*(2*pi./ld)/ctho),-1);
n=retsqrt(n_ext^2+A./((2*pi./ld0).^2-(2*pi./ld).^2-1i*(2*pi./ld)/ctho),-1);% Modif 2 2013  bonne formule ?

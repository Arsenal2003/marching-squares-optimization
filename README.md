# marching-squares-optimization

Fiecare thread realizeaza secvential o parte din fiecare pas al algoritmului, calculeaza imaginea rescalata, creeaza grid-ul si asigneaza forma de contur specifica(rescale_image,sample_grid,march).
Pentru a impiedica citirea incorecta la un pas al algoritmului, am utilizat o bariera intre pasul de rescale_image si cel de sample_grid, cat si intre sample_grid si march pentru a asigura ca matricile grid si scaled_image au fost calculate in totalitate.    
 

threads_aux_info struct : 
  id -> id_ul fiecarui thread 
  NR_THREADS -> nr de threaduri care paralelizeaza alg
  image -> pointer catre imaginea primita ca parametru 
  scaled_image -> pointer catre imaginea rescalata, daca imaginea initiala a fost prea mare 
  contour_map -> pointer catre countour map 
  grid -> pointer catre grid   
  bariera -> pointer catre bariera folosita de threaduri






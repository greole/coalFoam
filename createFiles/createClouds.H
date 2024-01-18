clock_t t_begin_pc1 = std::clock();
Info<< "\nConstructing coal cloud1" << endl;

coalCloud coalParcels (
    "coalCloud1",
    rho,
    U,
    g,
    slgThermo
);

probesContainer particle_probes = probesContainer(runTime, mesh, coalParcels);

clock_t t_end_pc1 = std::clock();
double elapsed_secs_pc1 = double(t_end_pc1 - t_begin_pc1) / CLOCKS_PER_SEC;
Info <<  "===  [Done]: " << elapsed_secs_pc1 << " s ===\n" << endl;

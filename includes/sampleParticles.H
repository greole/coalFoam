clock_t t_begin_ps = std::clock();
Info <<  "=== Start Particle Sampling ===" << endl;

particle_probes.sample(coalParcels);

clock_t t_end_ps = std::clock();
double elapsed_secs_ps = double(t_end_ps - t_begin_ps) / CLOCKS_PER_SEC;
Info <<  "=== Particle Sampling [Done]: " << elapsed_secs_ps << " s ===\n" << endl;

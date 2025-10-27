# First, check which packages renv thinks are installed
installed_packages <- renv::dependencies()$Package
unavailable <- c("AssocTests", "tess3r", "LEA")  # Add other problematic packages here

# Snapshot excluding problematic packages
renv::snapshot(exclude = unavailable)

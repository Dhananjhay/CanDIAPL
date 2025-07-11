import os
import pyvo

# TAP service endpoint
service = pyvo.dal.TAPService("https://gea.esac.esa.int/tap-server/tap")

query = """
SELECT TOP 100
  source_id, ra, dec, parallax, phot_g_mean_mag,
  DISTANCE(36.033303, 18.25795, ra, dec) AS ang_sep
FROM gaiadr3.gaia_source
WHERE DISTANCE(36.033303, 18.25795, ra, dec) < 0.1
ORDER BY ang_sep ASC
"""

results = service.search(query)
table = results.to_table()
output_file = snakemake.output["out_path"]

# Ensure output directory exists
os.makedirs(os.path.dirname(output_file), exist_ok=True)

# Write to CSV
table.write(output_file, format="ascii.csv", overwrite=True)
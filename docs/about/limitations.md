# Limitations and known issues

- **Very large molecules are skipped in some workflows.** A molecule with more
  than 30 non-ring single bonds is not mutated, and a molecule with more than
  100 hydrogen atoms is not processed by grow or link.

- **Database compatibility depends on the RDKit version.** Context
  canonicalization uses RDKit's SMILES representation. A change in how RDKit
  writes SMILES can make a database built with one RDKit version unusable from a
  different version. We did not observed this in the past, but recommend to pin 
  RDKit (`>=2025.3.5`) when sharing databases across environments.

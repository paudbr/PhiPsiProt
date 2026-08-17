# Imagen de métricas estructurales: PyRosetta + propka + biopython + freesasa + tmtools
#
# Construir:
#   docker build -t phipsiprot/pyrosetta:dev-metrics -f Dockerfile.metrics .

FROM phipsiprot/pyrosetta:dev

# Python 3.12 marca el entorno como "externally managed" (PEP 668),
# por eso --break-system-packages. Instala en el Python del sistema (donde vive pyrosetta).
#   propka    -> carga neta por pKa de entorno
#   biopython -> parseo PDB/CIF, pLDDT, cisteínas, secuencia para TM-score
#   freesasa  -> SASA para el SAP score
#   tmtools   -> TM-score real (algoritmo de Zhang)
RUN pip install --no-cache-dir --break-system-packages \
        propka \
        biopython \
        freesasa \
        tmtools

# Verificación en build (falla si algo no importa EN EL MISMO python)
RUN python3 -c "import pyrosetta, propka, Bio, freesasa, tmtools, numpy; print('all metric deps OK')"
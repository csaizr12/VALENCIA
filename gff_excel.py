import pandas as pd
import re
import os

# Rutas
input_file = 'test_valencia_araport11/Athaliana_447_Araport11.gene_exons_with_evidence_features.gff3'
output_file = 'mRNAs_perfectos_con_exones.xlsx'

def extraer(attr, clave):
    m = re.search(fr'{clave}=([^;]+)', attr)
    return m.group(1) if m else None

print(f"Analizando archivo: {input_file}...")

mRNAs = []
conteo_exones = {}

with open(input_file, 'r') as f:
    for line in f:
        if line.startswith('#') or not line.strip():
            continue
            
        cols = line.split('\t')
        if len(cols) < 9:
            continue
        
        tipo = cols[2]
        attrs = cols[8]
        
        if tipo == 'mRNA':
            id_mrna = extraer(attrs, 'ID')
            dist = extraer(attrs, 'transcripts_edit_distance')
            
            # Solo guardamos los de distancia 0
            if dist == '0':
                mRNAs.append({
                    'Cromosoma': cols[0],
                    'ID_mRNA': id_mrna,
                    'Class_Code': extraer(attrs, 'transcripts_class_code'),
                    'Edit_Distance': 0,
                    'Inicio': cols[3],
                    'Fin': cols[4]
                })
        
        elif tipo == 'exon':
            parent = extraer(attrs, 'Parent')
            if parent:
                # Contamos cuántos exones tiene cada Parent ID
                conteo_exones[parent] = conteo_exones.get(parent, 0) + 1

# Unimos la información: mRNA + su conteo de exones
for mrna in mRNAs:
    id_actual = mrna['ID_mRNA']
    mrna['Num_Exones'] = conteo_exones.get(id_actual, 0)

# Crear DataFrame y guardar
df = pd.DataFrame(mRNAs)

if not df.empty:
    df.to_excel(output_file, index=False)
    print(f"\n✅ Proceso completado.")
    print(f"Se han encontrado {len(df)} mRNAs con distancia 0.")
    print(f"Archivo generado: {output_file}")
    
    # Resumen rápido en terminal
    print("\n--- Resumen por Class Code ---")
    print(df.groupby('Class_Code')['Num_Exones'].agg(['count', 'sum']).rename(columns={'count': 'Num_mRNAs', 'sum': 'Total_Exones'}))
else:
    print("⚠️ No se encontró ningún mRNA con distancia de editado 0.")
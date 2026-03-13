import pandas as pd
import re

def extraer_valor(atributo_str, clave):
    """Busca una clave exacta y devuelve su valor."""
    match = re.search(fr'{clave}=([^;]+)', atributo_str)
    return match.group(1) if match else "N/A"

def convertir_gff():
    input_file = 'Athaliana_447_Araport11.gene_exons_with_evidence_features.gff3' 
    output_file = 'revision_anotacion.xlsx'
    
    registros = []

    print(f"Leyendo {input_file}...")

    with open(input_file, 'r') as f:
        for linea in f:
            # Saltamos encabezados
            if linea.startswith('#') or not linea.strip():
                continue
            
            campos = linea.strip().split('\t')
            if len(campos) < 9:
                continue
            
            tipo_feature = campos[2]
            
            # Filtramos solo lo que te interesa
            if tipo_feature in ['mRNA', 'exon']:
                atributos = campos[8]
                
                # Extraemos los datos clave
                edit_dist = extraer_valor(atributos, 'transcripts_edit_distance')
                class_code = extraer_valor(atributos, 'transcripts_class_code')
                id_name = extraer_valor(atributos, 'ID')
                parent = extraer_valor(atributos, 'Parent')

                registros.append({
                    'Cromosoma': campos[0],
                    'Tipo': tipo_feature,
                    'Inicio': int(campos[3]),
                    'Fin': int(campos[4]),
                    'ID': id_name,
                    'Parent': parent,
                    'Class_Code': class_code,
                    'Edit_Distance': edit_dist
                })

    # Creamos el DataFrame
    df = pd.DataFrame(registros)

    # Convertimos Edit_Distance a número (los N/A serán NaN)
    df['Edit_Distance'] = pd.to_numeric(df['Edit_Distance'], errors='coerce')

    # Guardar a Excel
    df.to_excel(output_file, index=False)
    
    # --- Pequeño reporte por terminal ---
    print("\n" + "="*40)
    print("¡CONVERSIÓN COMPLETADA!")
    print(f"Total de registros: {len(df)}")
    print(f"mRNAs: {len(df[df['Tipo'] == 'mRNA'])}")
    print(f"Exones: {len(df[df['Tipo'] == 'exon'])}")
    print(f"Registros con Edit Distance = 0: {len(df[df['Edit_Distance'] == 0])}")
    print("="*40)
    print(f"Archivo generado: {output_file}")

if __name__ == "__main__":
    convertir_gff()
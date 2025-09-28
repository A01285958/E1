#Lectura de Archivos 
"""
Lee un archivo de secuencia y solo devuelve las
secuencias de ADN (A, C, G, T, N). Ignora encabezados
y otros caracteres
 """
def leer_secuencia(ruta):
    secuencia = ""
    with open(ruta, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip().upper()
            #Solo secuencias validas que contengan ACGTN
            if all(c in "ACGTN" for c in line): 
                secuencia += line
    return secuencia

"""
Se calcula el reverso complementario de la secuencia 
para asegurar detectar la posicion del gen aunque este escrito en la 
otra hebra del ADN, el ADN es doble cadena
"""
def reverso_complementario(seq):
    #Diccionario de bases complementarias
    complementos = {
        "A" : "T",
        "T" : "A",
        "C" : "G",
        "G" : "C",
        "N" : "N", #N se queda igual porque es base desconocida
    }

    #Se genera la cadena complementaria
    comp = ""
    for base in seq:
        comp += complementos.get(base, "N") # si aparece algo raro, lo convierte en N
    
    #Regresa la cadena al reves
    return comp[::-1]

"""
Busca todas las apariciones(1-based como en articulos cientficos) de un gen dentro de la
hebra. Devuelve lista de tuplas (inicio, fin)
"""

def buscar_todas(hebra, gen):
    indices = []
    n = len(hebra)
    m = len(gen)
    for i in range(n - m + 1):
        # compara el fragmento
        if hebra[i:i + m] == gen:
            # Se suma +1 al inicio para pasarlo a 1-based
            indices.append((i+1, i+m)) #Fin inclusivo
    return indices

def resultados_indices(nombre_genoma, genoma, genes):
    print(f"\nIndices de aparicion en {nombre_genoma}")
    for nombre, sec in genes.items():
        directo = buscar_todas(genoma, sec)
        complemento = buscar_todas(genoma, reverso_complementario(sec))
        print(f"{nombre}: primeros 12 = {sec[:12]}")
        print(f"Directo {directo if directo else "ninguno"}")
        print(f"Reverso Complementario {complemento if complemento else "ninguno"}")

def main():
    genomas = {
        "Wuhan_2019": "SARS-COV-2-MN908947.3.txt",
        "Texas_2020":"SARS-COV-2-MT106054.1.txt"
    }
    genes = {
        "Gen M": "gen-M.txt",
        "Gen ORF1ab": "gen-ORF1AB.txt",
        "Gen S": "gen-S.txt"
    }

    #Leer genomas
    secuencias_genomas = {n: leer_secuencia(r) for n, r in genomas.items()}
    #Leer genes
    secuencias_genes = {n: leer_secuencia(r) for n, r in genes.items()}

    #Reporte para Wuhan
    resultados_indices("Wuhan_2019", secuencias_genomas["Wuhan_2019"], secuencias_genes)

if __name__ == "__main__":
    main()

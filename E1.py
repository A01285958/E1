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

def LPS(gen):
  lps = [0] * len(gen)
  j = 0 #Longitud del prefijo actual
  i = 1
  while i < len(gen):
    if gen[i] == gen[j]:
      j += 1
      lps[i] = j
      i += 1
    else:
      if j != 0:
        j = lps[j -1] #Retrocede al mejor prefijo anterior
      else:
        lps[i] = 0
        i += 1
  return lps

def buscar_todas(secuencia, gen):
  if not gen:
    return []
  pos = []
  n = len(secuencia)
  m = len(gen)
  lps = LPS(gen)
  i = 0 #Indice en secuencia
  j = 0 #Indice en gen
  while (i < n):
    #Si coinciden, ambos punteros avanzan
    if(secuencia[i] == gen[j]):
      i += 1
      j += 1
    #Si no hay coincidencia j retrocede
    else:
      if(j > 0):
        j = lps[j - 1] #usa LPS para saltar comparaciones
      else:
        i += 1

    if(j == m):
      inicio = i - m
      fin = i - 1
      pos.append((inicio, fin))  #Coincidencia completa
      j = lps[j - 1]
  return pos


def resultados_indices(nombre_genoma, genoma, genes):
    print(f"\nIndices de aparicion en {nombre_genoma}")
    for nombre, sec in genes.items():
        directo = buscar_todas(genoma, sec)
        complemento = buscar_todas(genoma, reverso_complementario(sec))
        print(f"{nombre}: primeros 12 = {sec[:12]}")
        print(f"Directo {directo if directo else "ninguno"}")
        print(f"Reverso Complementario {complemento if complemento else "ninguno"}")

"""
Encuentra el palindromo mas largo en una secuencia
Utilizando la funcion transform para convertir los 3 genes
Utilizando el algoritmo de Manacher para encontrar el palindromo
"""
def transform(s: str):
    # Devolvemos lista (tu código ya indexa como lista)
    texto = ['@', '$']
    for c in s:
        texto.append(c)
        texto.append('$')
    texto.append('#')
    return texto

def Manacher(gen):
    texto = transform(gen)
    e = len(texto)
    P = [0] * e
    centro, limite = 0, 0

    for i in range(1, e - 1):
        if i < limite:
            sim = 2 * centro - i
            P[i] = min(limite - i, P[sim])

        # expandir alrededor de i
        gap = P[i] + 1
        # Con sentinelas, esta comparación se detiene antes de salirse
        while texto[i - gap] == texto[i + gap]:
            P[i] += 1
            gap += 1

        # actualizar centro y límite
        if i + P[i] > limite:
            limite = i + P[i]
            centro = i

    # localizar mejor centro (ojo: NO reutilizar la i anterior)
    best_i = 1
    best_len = 0
    for k in range(1, e - 1):
        if P[k] > best_len:
            best_len = P[k]
            best_i = k

    # índices en la cadena original
    start = (best_i - best_len) // 2
    end = start + best_len - 1
    return start, end

def resultados_palindromos(genes: dict):
    orden = ["Gen M", "Gen S", "Gen ORF1ab"]
    for nombre in orden:
        sec = genes[nombre]
        i, j = Manacher(sec)
        pal = sec[i:j+1] if j >= i >= 0 else ""
        print(f"{nombre}: longitud={len(pal)}, palindromo='{pal}'")

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
    resultados_palindromos(secuencias_genes)
    

if __name__ == "__main__":
    main()

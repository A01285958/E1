# Lectura de Archivos
"""
Lee un archivo de secuencia y solo devuelve las
secuencias de ADN (A, C, G, T, N). Ignora encabezados
y otros caracteres
"""


def leer_secuencia(ruta):
    secuencia = ""
    with open(ruta, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip().upper()
            # Solo secuencias validas que contengan ACGTN
            if all(c in "ACGTN" for c in line):
                secuencia += line
    return secuencia


"""
Lee un archivo de proteinas en formato FASTA
Cada linea que empieza con > es el nombre de la proteina
y la siguiente linea contiene la secuencia
Devuelve un diccionario: {nombre: secuencia}
"""


def leer_proteinas(ruta):
    proteinas = {}
    with open(ruta, "r", encoding="utf-8", errors="ignore") as f:
        lineas = [linea.strip() for linea in f if linea.strip()]
    for i in range(len(lineas)):
        if lineas[i].startswith(">"):
            nombre = lineas[i][1:]  # Quita el >
            if i + 1 < len(lineas):
                secuencia = lineas[i + 1]
                proteinas[nombre] = secuencia
    return proteinas


"""
Se calcula el reverso complementario de la secuencia 
para asegurar detectar la posicion del gen aunque este escrito en la 
otra hebra del ADN, el ADN es doble cadena
"""


def reverso_complementario(seq):
    # Diccionario de bases complementarias
    complementos = {
        "A": "T",
        "T": "A",
        "C": "G",
        "G": "C",
        "N": "N",  # N se queda igual porque es base desconocida
    }

    # Se genera la cadena complementaria
    comp = ""
    for base in seq:
        comp += complementos.get(base, "N")  # si aparece algo raro, lo convierte en N

    # Regresa la cadena al reves
    return comp[::-1]


"""
Busca todas las apariciones de un gen dentro de la
hebra. Devuelve lista de tuplas (inicio, fin)
"""


def LPS(gen):
    lps = [0] * len(gen)
    j = 0  # Longitud del prefijo actual
    i = 1
    while i < len(gen):
        if gen[i] == gen[j]:
            j += 1
            lps[i] = j
            i += 1
        else:
            if j != 0:
                j = lps[j - 1]  # Retrocede al mejor prefijo anterior
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
    i = 0  # Indice en secuencia
    j = 0  # Indice en gen
    while i < n:
        # Si coinciden, ambos punteros avanzan
        if secuencia[i] == gen[j]:
            i += 1
            j += 1
        # Si no hay coincidencia j retrocede
        else:
            if j > 0:
                j = lps[j - 1]  # usa LPS para saltar comparaciones
            else:
                i += 1

        if j == m:
            inicio = i - m
            fin = i - 1
            pos.append((inicio, fin))  # Coincidencia completa
            j = lps[j - 1]
    return pos


def resultados_indices(nombre_genoma, genoma, genes):
    print(f"\nIndices de aparicion en {nombre_genoma}")
    for nombre, sec in genes.items():
        directo = buscar_todas(genoma, sec)
        complemento = buscar_todas(genoma, reverso_complementario(sec))

        print(f"{nombre}: primeros 12 = {sec[:12]}")
        print(f"Directo {directo if directo else "ninguno"}")
        if complemento:
            print(f"Reverso Complementario {complemento}")


# ==== Proteinas ====
# ====================
# Tabla genética (DNA)
# ====================
CODON_TABLE = {
    "TTT": "F",
    "TTC": "F",
    "TTA": "L",
    "TTG": "L",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "ATT": "I",
    "ATC": "I",
    "ATA": "I",
    "ATG": "M",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "AGT": "S",
    "AGC": "S",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "TAT": "Y",
    "TAC": "Y",
    "CAT": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "AAT": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "GAT": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "TGT": "C",
    "TGC": "C",
    "TGG": "W",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "AGA": "R",
    "AGG": "R",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
    "TAA": "*",
    "TAG": "*",
    "TGA": "*",
}


def traducir_codon(codon):
    # Diccionario codon -> aminoácido (tabla estándar)
    return CODON_TABLE.get(codon, "X")


# Traduce en un marco (offset 0/1/2). N→X si el codón no está en la tabla
def traducir_hebra(seq, offset=0):
    aa = []
    for i in range(offset, len(seq) - 2, 3):
        codon = seq[i : i + 3]
        aa.append(CODON_TABLE.get(codon, "X"))
    return "".join(aa)


# ==== Traduccion en 6 marcos =====
def traducir_6marcos(genoma):
    """Devuelve dict con genoma traducido a aminoacidos"""
    out = {}
    # Hebra directa
    for f in range(3):
        out[f"+frame{f+1}"] = traducir_hebra(genoma, f)
    # Hebra inversa
    rev = reverso_complementario(genoma)
    for f in range(3):
        out[f"-frame{f+1}"] = traducir_hebra(rev, f)
    return out


# ==== Mapeo de posiciones de aminoacidos a coordendas NT ====
def map_AA_NT(gen_len, frame_key, aa_i, aa_f):
    """Devuelve indices de incion en la secuencia original"""
    plus = frame_key.startswith("+")
    fnum = int(frame_key[-1]) - 1  # offset 0,1,2
    aa_len = aa_f - aa_i + 1

    # Para casos especiales (como la proteina 11)
    # restar 1 solo si NO es +frame1
    if plus:
        nt_start = fnum + (aa_i * 3)
        nt_end = nt_start + aa_len * 3 - 1
        return (nt_start, nt_end)
    else:
        # En la hebra reversed complementario el indice 0 corresponde a nt og gen_len-1
        start_rc = fnum + (aa_i * 3)
        end_rc = start_rc + aa_len * 3 - 1
        nt_start = (gen_len - 1) - end_rc
        nt_end = (gen_len - 1) - start_rc
        return (nt_start, nt_end)


def traducir_desde(genoma, frame_key, nt_start, aa_len):
    # Traduce aa_len AA desde nt_start en el mismo marco
    frag = genoma[nt_start : nt_start + aa_len * 3]
    if frame_key.startswith("+"):
        return traducir_hebra(frag, 0)
    else:
        frag_rc = reverso_complementario(frag)
        return traducir_hebra(frag_rc, 0)


# === Busqueda de Proteinas ====
# nt_end = nt_start + len(AA)*3 - 1 (un final estimado, puede cruzar frameshift)
def buscar_proteinas(genoma, proteinas_dict, k=4, verif_len=120):
    # Diccionario donde se almacenaran los resultaods:{nombre_prot:[(inicio, fin, primeros4, codones)]}
    resultados = {}

    # Traducción del genoma completo en los 6 marcos de lectura
    # +frame1, +frame2, +frame3 - lectura directa
    # -frame1, -frame2, -frame3 - reverso complementario
    marcos = traducir_6marcos(genoma)
    n = len(genoma)  # Longitud del genoma (virus)

    # Itera por cada proteina del seq-proteins.txt
    for nombre, aa in proteinas_dict.items():
        if not aa:
            continue  # Si esta vacia, se salta

        # Se toma un prefijo (primeros k aminoacidos) como patron para busqueda rapida
        aa_pref = aa[: k + 1]
        busquedas = []

        # Se recorre cada uno de los 6 marcos de lectura traducidos
        for frame_key, aa_seq in marcos.items():
            # Busca todas las aparicion de aa_pref en la secuencia traducida
            matches = buscar_todas(aa_seq, aa_pref)

            # Por cada coincidencia se obtiene una posicion dentro del marco
            for aa_i, aa_f in matches:  # Posciones en aa

                # Convetir las posiciones de aa a coordenadas de nucleotidos
                nt_start, _ = map_AA_NT(n, frame_key, aa_i, aa_f)

                # Verificacion: comparar una extension mas larga
                # Traduce un fragmento las largo
                L = min(verif_len, len(aa))
                aa_ext = traducir_desde(genoma, frame_key, nt_start, L)
                # Compara caracter por caracter hasta que encuentra una diferencia
                match_len = 0
                for x, y in zip(aa_ext, aa[:L]):
                    if x != y:
                        break
                    match_len += 1
                # Guarda coincidencia (longitud de coincidencia, marco, posición de inicio)
                busquedas.append((match_len, frame_key, nt_start))

        # Si se encontraron coincidencias
        if busquedas:
            # Ordenar para priorizar el match mas largo
            busquedas.sort(reverse=True)  # Mayor match_len primero
            mejor = busquedas[0]
            match_len, frame_key, nt_start = mejor
            # Codones asociados (primeros 12)
            codones = genoma[nt_start : nt_start + 12]
            if not frame_key.startswith("+"):
                codones = reverso_complementario(codones)
            # Fin estimado si toda la proteina estuviera en un marco
            nt_end = nt_start + len(aa) * 3 - 1
            resultados[nombre] = [(nt_start, nt_end, aa[:4], codones)]
    return resultados


def resultados_proteinas(resultados):
    for nombre, lista in resultados.items():
        for nt_start, nt_end, first4, codones in lista:
            print(f"Proteina: {nombre}")
            print(f" Genoma: {nt_start} - {nt_end}")
            print(f" 1eros 4 AA: {first4}")
            print(f" Codones: {codones}")
            print()


"""
Encuentra el palindromo mas largo en una secuencia
Utilizando la funcion transform para convertir los 3 genes
Utilizando el algoritmo de Manacher para encontrar el palindromo
"""


def transform(s: str):
    # Devolvemos lista (tu código ya indexa como lista)
    texto = ["@", "$"]
    for c in s:
        texto.append(c)
        texto.append("$")
    texto.append("#")
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
        pal = sec[i : j + 1] if j >= i >= 0 else ""
        print(f"{nombre}: longitud={len(pal)}, palindromo='{pal}'")


def main():
    genomas = {
        "Wuhan_2019": "SARS-COV-2-MN908947.3.txt",
        "Texas_2020": "SARS-COV-2-MT106054.1.txt",
    }
    genes = {"Gen M": "gen-M.txt", "Gen ORF1ab": "gen-ORF1AB.txt", "Gen S": "gen-S.txt"}

    # Leer genomas
    secuencias_genomas = {n: leer_secuencia(r) for n, r in genomas.items()}
    # Leer genes
    secuencias_genes = {n: leer_secuencia(r) for n, r in genes.items()}

    proteinas = leer_proteinas("seq-proteins.txt")

    # Reporte para Wuhan
    print("Genes en Virus\n")
    resultados_indices("Wuhan_2019", secuencias_genomas["Wuhan_2019"], secuencias_genes)
    print()
    print("Palindromos mas largos\n")
    resultados_palindromos(secuencias_genes)
    print()
    print("Proteinas\n")
    res = buscar_proteinas(secuencias_genomas["Wuhan_2019"], proteinas)
    resultados_proteinas(res)

    # for nombre, sec in list(proteinas.items())[:3]:
    #     print(nombre, sec[:30], "...")


if __name__ == "__main__":
    main()

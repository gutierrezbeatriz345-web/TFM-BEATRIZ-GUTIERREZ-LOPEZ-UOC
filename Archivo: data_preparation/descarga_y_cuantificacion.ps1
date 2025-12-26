# ===== CONFIGURACIÓN =====

# Carpeta donde guardar los .sra
$dest = "D:\"

# Lista de SRR a descargar
$srrs = @(
  "SRR13435062",
  "SRR13435063",
  "SRR13435064",
  "SRR13435058"
)

# Crea carpeta si no existe
New-Item -ItemType Directory -Force -Path $dest | Out-Null

# ===== DESCARGA =====

foreach ($srr in $srrs) {
    Write-Host ">> Descargando $srr ..."
    prefetch $srr --output-directory $dest
    Write-Host " $srr descargado en $dest"
}

Write-Host " Descargas completadas."

# ===== CONVERSIÓN a FASTQ =====

# Carpeta donde están los .sra (descargados en D:)
$src = "D:\"

# Carpeta de salida para los FASTQ
$dest = "C:\fastq_out"

# Carpeta temporal (necesaria para fasterq-dump)
$tmp = "C:\tmp_sra"

# Crea las carpetas si no existen
New-Item -ItemType Directory -Force -Path $dest | Out-Null
New-Item -ItemType Directory -Force -Path $tmp  | Out-Null

# Lista de SRR a convertir (los que ya tienes descargados)
$srrs = @(
  "SRR13435063",
  "SRR13435073",
  "SRR13435074",
  "SRR13435075"
)

# Conversión SRA -> FASTQ
foreach ($srr in $srrs) {
    Write-Host ">> Convirtiendo $srr ..."
    $sraFile = Join-Path $src "$srr\$srr.sra"

    if (Test-Path $sraFile) {
        fasterq-dump $sraFile -O $dest -t $tmp -e 4 -p
        Write-Host " $srr convertido correctamente.`n"
    } else {
        Write-Warning " No se encontró el archivo $sraFile"
    }
}

Write-Host " Conversión finalizada. Los FASTQ están en $dest"

# ===== CUANTIFICACIÓN con Salmon =====

# Montar la unidad E: en WSL (Ubuntu)
wsl sudo mkdir -p /mnt/e
wsl sudo mount -t drvfs E: /mnt/e

# Ejecutar Salmon para cada muestra SRR (desde WSL a través de PowerShell)
foreach ($srr in $srrs) {
    Write-Host ">> Cuantificando $srr con Salmon..."
    # Ejecuta Salmon en WSL (índice pre-construido en /mnt/e/ref/salmon_index)
    wsl salmon quant -i /mnt/e/ref/salmon_index -l IU `
        -1 "/mnt/c/fastq_out/${srr}_1.fastq" `
        -2 "/mnt/c/fastq_out/${srr}_2.fastq" `
        -p 4 -o "/mnt/e/quant_${srr}"
    Write-Host " $srr cuantificado correctamente."
}

Write-Host " Cuantificación completada."

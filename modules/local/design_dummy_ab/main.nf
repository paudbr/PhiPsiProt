process DESIGN_DUMMY {
    tag "$mode"
    container 'rfantibody_local:1.0.0'

    // Descomenta la siguiente línea si tu clúster usa Slurm y requiere una etiqueta para GPUs
    // label 'process_gpu'

    input:
    val mode

    output:
    path "resultado_test_env.txt", emit: report
    stdout emit: stdout_log

    script:
    """
    echo "=================================================" > resultado_test_env.txt
    echo "  EJECUTANDO PROCESO DUMMY DE VALIDACIÓN" >> resultado_test_env.txt
    echo "=================================================" >> resultado_test_env.txt
    echo "Modo de pipeline solicitado: $mode" >> resultado_test_env.txt
    echo "Fecha de ejecución: \$(date)" >> resultado_test_env.txt
    echo "Usuario dentro del contenedor: \$(whoami)" >> resultado_test_env.txt
    echo "" >> resultado_test_env.txt

    echo "--- [1/3] VERIFICANDO PATH Y ENTORNO ---" >> resultado_test_env.txt
    echo "Ubicación de Python: \$(which python3)" >> resultado_test_env.txt
    echo "Versión de Python: \$(python3 --version)" >> resultado_test_env.txt
    echo "" >> resultado_test_env.txt

    echo "--- [2/3] PRUEBA DE LIBRERÍAS CRÍTICAS ---" >> resultado_test_env.txt
    
    python3 -c "
import sys, torch
print('-> Ruta real del intérprete:', sys.executable)
print('-> Versión de PyTorch:', torch.__version__)
print('-> ¿CUDA (GPU) disponible?:', torch.cuda.is_available())
if torch.cuda.is_available():
    print('   -> Nombre de la GPU:', torch.cuda.get_device_name(0))
else:
    print('   ⚠️ ALERTA: PyTorch no detecta la GPU dentro de Nextflow.')

try:
    import dgl
    print('-> Versión de DGL:', dgl.__version__)
    print('\n✔ EXCELENTE: Todo el entorno virtual de RFantibody está intacto.')
except Exception as e:
    print('\n❌ ERROR CRÍTICO: No se pudo importar DGL:', e)
    sys.exit(1)
" >> resultado_test_env.txt 2>&1

    echo "" >> resultado_test_env.txt
    echo "--- [3/3] REPOSITORIO INTERNO ---" >> resultado_test_env.txt
    echo "Archivos en /opt/rfantibody:" >> resultado_test_env.txt
    ls -l /opt/rfantibody | head -n 5 >> resultado_test_env.txt

    # Mostrar el reporte final directamente en la consola de Nextflow
    cat resultado_test_env.txt
    """
}
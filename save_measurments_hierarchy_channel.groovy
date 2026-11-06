import org.locationtech.jts.geom.Coordinate
import org.locationtech.jts.geom.prep.PreparedGeometryFactory

// ============================================================
// CONFIGURACIÓN DE LAS CLASES DE ANOTACIÓN
// ============================================================

/*
 * annotationClassName:
 * Nombre real de la clase de la anotación en QuPath.
 *
 * detectionClassName:
 * Nombre que se asignará a las detecciones y que aparecerá
 * en la columna Classification del archivo exportado.
 *
 * IMPORTANTE:
 * El orden determina la prioridad.
 */
def annotationClassPriority = [
    [
        annotationClassName: 'Lesion_dermis',
        detectionClassName : 'Annotation (Lesion_dermis)'
    ],
    [
        annotationClassName: 'Lesion_epidermis',
        detectionClassName : 'Annotation (Lesion_epidermis)'
    ],
    [
        annotationClassName: 'Lesion',
        detectionClassName : 'Annotation (Lesion)'
    ],
    [
        annotationClassName: 'dermis',
        detectionClassName : 'Annotation (dermis)'
    ],
    [
        annotationClassName: 'epidermis',
        detectionClassName : 'Annotation (epidermis)'
    ]
]

// ============================================================
// CONFIGURACIÓN DE LA CLASIFICACIÓN
// ============================================================

/*
 * true:
 * Las detecciones fuera de todas las anotaciones configuradas
 * quedarán sin clasificación.
 *
 * false:
 * Conservarán la clasificación que tuvieran anteriormente.
 */
boolean resetDetectionsOutsideAnnotations = true

/*
 * true:
 * Si la detección es una célula y tiene núcleo, se utiliza
 * el centroide del núcleo.
 *
 * false:
 * Se utiliza el centroide de la célula o detección completa.
 */
boolean useNucleusCentroid = true

/*
 * Añadir una medición numérica adicional:
 *
 * 0 = Fuera de las regiones
 * 1 = Annotation (Lesion_dermis)
 * 2 = Annotation (Lesion_epidermis)
 * 3 = Annotation (Lesion)
 * 4 = Annotation (dermis)
 * 5 = Annotation (epidermis)
 */
boolean addRegionCodeMeasurement = true

String regionMeasurementName = 'Annotation region code'

// Mostrar progreso cada 5.000 detecciones
int progressInterval = 5000

// ============================================================
// CONFIGURACIÓN DE LA EXPORTACIÓN
// ============================================================

// Carpeta de salida dentro del proyecto
String outputFolderName = 'measurements_cellpose_process'

// Nombre final del archivo
String outputFileSuffix = '_single_cell.csv'

// Nombre temporal utilizado durante la exportación
String temporaryFileSuffix =
    '_single_cell_original_channel_names.csv'

/*
 * true:
 * Conserva también el archivo con los nombres AF originales.
 *
 * false:
 * Elimina el archivo temporal después de crear el CSV final.
 */
boolean keepOriginalCSV = false

// ============================================================
// CORRESPONDENCIA ENTRE CANALES ANTIGUOS Y NUEVOS
// ============================================================

/*
 * Esta correspondencia cambia únicamente los encabezados
 * del archivo exportado.
 *
 * Los valores numéricos de las mediciones no se modifican.
 */
def channelNameMapping = [
    'AF405': 'DAPI',
    'AF488': 'K14',
    'AF555': 'MDK',
    'AF647': 'PTN',
    'AF750': 'S100B'
]

// ============================================================
// OBTENER IMAGEN, JERARQUÍA Y OBJETOS
// ============================================================

def imageData = getCurrentImageData()

if (imageData == null) {
    println 'ERROR: No hay ninguna imagen abierta.'
    return
}

def hierarchy = imageData.getHierarchy()
def server = imageData.getServer()

def detections = getDetectionObjects()
def annotations = getAnnotationObjects()

if (detections == null || detections.isEmpty()) {
    println 'ERROR: No se encontraron detecciones.'
    return
}

if (annotations == null || annotations.isEmpty()) {
    println 'ERROR: No se encontraron anotaciones.'
    return
}

// ============================================================
// OBTENER EL NOMBRE DE LA IMAGEN
// ============================================================

def projectEntry = getProjectEntry()

String imageName

if (projectEntry != null) {
    imageName = projectEntry.getImageName()
} else {
    imageName = server.getMetadata().getName()
}

/*
 * Eliminar extensiones conocidas y sustituir caracteres
 * problemáticos para crear un nombre de archivo seguro.
 */
String safeName = imageName
    .replaceFirst(
        '(?i)\\.(ome\\.tiff|ome\\.tif|tiff|tif|czi|lif)$',
        ''
    )
    .replaceAll('[^a-zA-Z0-9_.-]', '_')
    .replaceAll('_+', '_')
    .replaceAll('^_+', '')
    .replaceAll('_+$', '')

if (safeName == null || safeName.trim().isEmpty()) {
    safeName = 'image'
}

// ============================================================
// PREPARAR LAS RUTAS DE SALIDA
// ============================================================

def outputDir = buildFilePath(
    PROJECT_BASE_DIR,
    outputFolderName
)

mkdirs(outputDir)

// Archivo final con nombres DAPI, K14, MDK, PTN y S100B
def detectionMeasurementsPath = buildFilePath(
    outputDir,
    safeName + outputFileSuffix
)

// Archivo temporal con los nombres AF originales
def temporaryMeasurementsPath = buildFilePath(
    outputDir,
    safeName + temporaryFileSuffix
)

// Evitar que el archivo final y el temporal sean el mismo
if (detectionMeasurementsPath == temporaryMeasurementsPath) {
    println 'ERROR: La ruta temporal coincide con la ruta final.'
    return
}

// ============================================================
// FUNCIÓN AUXILIAR PARA OBTENER EL NOMBRE DE CLASE
// ============================================================

def getObjectClassName = { pathObject ->

    def pathClass = pathObject.getPathClass()

    if (pathClass == null) {
        return null
    }

    return pathClass.toString()
}

// ============================================================
// SELECCIONAR ANOTACIONES RELEVANTES
// ============================================================

def relevantAnnotations = annotations.findAll { annotation ->

    String currentClassName =
        getObjectClassName(annotation)

    if (currentClassName == null) {
        return false
    }

    return annotationClassPriority.any { classDefinition ->

        String expectedClassName =
            classDefinition.annotationClassName as String

        currentClassName.equalsIgnoreCase(
            expectedClassName
        )
    }
}

// ============================================================
// COMPROBAR LAS CLASES ENCONTRADAS
// ============================================================

if (relevantAnnotations.isEmpty()) {

    println ''
    println 'ERROR: No se encontraron anotaciones con las clases:'

    annotationClassPriority.each { classDefinition ->

        println(
            '  Annotation (' +
            classDefinition.annotationClassName +
            ')'
        )
    }

    println ''
    println 'Clases reales encontradas en la imagen:'

    def availableClasses = annotations
        .collect { annotation ->
            getObjectClassName(annotation)
        }
        .findAll { className ->
            className != null
        }
        .unique()
        .sort()

    availableClasses.each { className ->
        println '  ' + className
    }

    return
}

// ============================================================
// PREPARAR GEOMETRÍAS DE LAS ANOTACIONES
// ============================================================

println ''
println '=========================================='
println 'PREPARACIÓN DE ANOTACIONES'
println '=========================================='
println 'Imagen: ' + imageName
println 'Detecciones: ' + detections.size()
println 'Anotaciones totales: ' + annotations.size()
println 'Anotaciones relevantes: ' + relevantAnnotations.size()
println ''

def preparedAnnotations = []

annotationClassPriority.eachWithIndex {
    classDefinition,
    priorityIndex ->

    String annotationClassName =
        classDefinition.annotationClassName as String

    String detectionClassName =
        classDefinition.detectionClassName as String

    def annotationsForClass =
        relevantAnnotations.findAll { annotation ->

            String currentClassName =
                getObjectClassName(annotation)

            return currentClassName != null &&
                currentClassName.equalsIgnoreCase(
                    annotationClassName
                )
        }

    println(
        'Annotation (' +
        annotationClassName +
        '): ' +
        annotationsForClass.size()
    )

    annotationsForClass.each { annotation ->

        def roi = annotation.getROI()

        if (roi == null || roi.isEmpty()) {

            println(
                'ADVERTENCIA: Se ignoró una Annotation (' +
                annotationClassName +
                ') con ROI vacío.'
            )

            return
        }

        try {

            def geometry = roi.getGeometry()

            if (geometry == null || geometry.isEmpty()) {
                return
            }

            /*
             * Preparar la geometría una sola vez acelera
             * las comprobaciones de contención.
             */
            def preparedGeometry =
                PreparedGeometryFactory.prepare(
                    geometry
                )

            preparedAnnotations.add([
                annotationClassName : annotationClassName,
                detectionClassName  : detectionClassName,

                /*
                 * Clase que recibirá la detección.
                 * Aparecerá literalmente en Classification.
                 */
                detectionPathClass  : getPathClass(
                    detectionClassName
                ),

                priority            : priorityIndex,
                regionCode          : priorityIndex + 1,

                minX                : roi.getBoundsX(),
                minY                : roi.getBoundsY(),

                maxX                : roi.getBoundsX() +
                                      roi.getBoundsWidth(),

                maxY                : roi.getBoundsY() +
                                      roi.getBoundsHeight(),

                z                   : roi.getZ(),
                t                   : roi.getT(),

                geometryFactory     : geometry.getFactory(),
                preparedGeometry    : preparedGeometry
            ])

        } catch (Exception exception) {

            println(
                'ADVERTENCIA: No se pudo preparar Annotation (' +
                annotationClassName +
                '). Mensaje: ' +
                exception.getMessage()
            )
        }
    }
}

if (preparedAnnotations.isEmpty()) {

    println ''
    println 'ERROR: No se pudo preparar ninguna geometría válida.'
    return
}

// Mantener explícitamente el orden de prioridad
preparedAnnotations.sort { first, second ->
    first.priority <=> second.priority
}

// ============================================================
// AGRUPAR ANOTACIONES POR PLANO Z/T
// ============================================================

/*
 * Cada detección solo se compara con anotaciones que están
 * en el mismo plano Z y T.
 */
def annotationsByPlane =
    preparedAnnotations.groupBy { annotationInfo ->

        annotationInfo.z + ':' + annotationInfo.t
    }

// ============================================================
// PREPARAR CONTADORES
// ============================================================

def classificationCounts = [:]

annotationClassPriority.each { classDefinition ->

    String detectionClassName =
        classDefinition.detectionClassName as String

    classificationCounts[detectionClassName] = 0
}

int outsideAnnotations = 0
int invalidDetections = 0
int changedClassifications = 0
int unchangedClassifications = 0
int containmentErrors = 0

def changedObjects = []

long startTime = System.currentTimeMillis()

// ============================================================
// CLASIFICAR LAS DETECCIONES
// ============================================================

println ''
println '=========================================='
println 'CLASIFICACIÓN DE DETECCIONES'
println '=========================================='
println 'Centroide nuclear: ' + useNucleusCentroid
println 'Reiniciar detecciones externas: ' +
    resetDetectionsOutsideAnnotations
println ''

detections.eachWithIndex {
    detection,
    detectionIndex ->

    def detectionROI = detection.getROI()

    /*
     * Para células se usa el núcleo si está disponible.
     */
    if (
        useNucleusCentroid &&
        detection.metaClass.respondsTo(
            detection,
            'getNucleusROI'
        )
    ) {

        def nucleusROI = detection.getNucleusROI()

        if (
            nucleusROI != null &&
            !nucleusROI.isEmpty()
        ) {
            detectionROI = nucleusROI
        }
    }

    // --------------------------------------------------------
    // DETECCIÓN SIN ROI VÁLIDO
    // --------------------------------------------------------

    if (
        detectionROI == null ||
        detectionROI.isEmpty()
    ) {

        invalidDetections++

        if (addRegionCodeMeasurement) {

            def measurementList =
                detection.getMeasurementList()

            measurementList.put(
                regionMeasurementName,
                0.0
            )

            measurementList.close()
        }

        return
    }

    // --------------------------------------------------------
    // OBTENER CENTROIDE Y PLANO
    // --------------------------------------------------------

    double centroidX = detectionROI.getCentroidX()
    double centroidY = detectionROI.getCentroidY()

    int detectionZ = detectionROI.getZ()
    int detectionT = detectionROI.getT()

    String planeKey =
        detectionZ + ':' + detectionT

    def candidateAnnotations =
        annotationsByPlane[planeKey]

    if (candidateAnnotations == null) {
        candidateAnnotations = []
    }

    def matchingAnnotation = null

    // --------------------------------------------------------
    // BUSCAR LA ANOTACIÓN QUE CONTIENE EL CENTROIDE
    // --------------------------------------------------------

    /*
     * La primera coincidencia gana según esta prioridad:
     *
     * 1. Lesion_dermis
     * 2. Lesion_epidermis
     * 3. Lesion
     * 4. dermis
     * 5. epidermis
     */
    for (annotationInfo in candidateAnnotations) {

        /*
         * Comprobación rápida usando el bounding box.
         */
        if (
            centroidX < annotationInfo.minX ||
            centroidX > annotationInfo.maxX ||
            centroidY < annotationInfo.minY ||
            centroidY > annotationInfo.maxY
        ) {
            continue
        }

        try {

            def centroidPoint =
                annotationInfo.geometryFactory.createPoint(
                    new Coordinate(
                        centroidX,
                        centroidY
                    )
                )

            /*
             * covers() incluye también los centroides situados
             * exactamente en el borde de la anotación.
             */
            if (
                annotationInfo.preparedGeometry.covers(
                    centroidPoint
                )
            ) {

                matchingAnnotation = annotationInfo
                break
            }

        } catch (Exception exception) {

            containmentErrors++

            /*
             * Para evitar miles de mensajes, solo se muestran
             * los primeros diez errores.
             */
            if (containmentErrors <= 10) {

                println(
                    'ADVERTENCIA: Error comprobando detección ' +
                    (detectionIndex + 1) +
                    ': ' +
                    exception.getMessage()
                )
            }
        }
    }

    // --------------------------------------------------------
    // DETERMINAR LA NUEVA CLASIFICACIÓN
    // --------------------------------------------------------

    def currentClass = detection.getPathClass()
    def newClass = null

    int regionCode = 0

    if (matchingAnnotation != null) {

        newClass =
            matchingAnnotation.detectionPathClass

        regionCode =
            matchingAnnotation.regionCode as int

        classificationCounts[
            matchingAnnotation.detectionClassName
        ]++

    } else {

        outsideAnnotations++

        if (!resetDetectionsOutsideAnnotations) {
            newClass = currentClass
        }
    }

    // --------------------------------------------------------
    // ASIGNAR LA CLASE A LA DETECCIÓN
    // --------------------------------------------------------

    if (currentClass != newClass) {

        detection.setPathClass(newClass)

        changedClassifications++
        changedObjects.add(detection)

    } else {

        unchangedClassifications++
    }

    // --------------------------------------------------------
    // AÑADIR EL CÓDIGO NUMÉRICO DE REGIÓN
    // --------------------------------------------------------

    if (addRegionCodeMeasurement) {

        def measurementList =
            detection.getMeasurementList()

        measurementList.put(
            regionMeasurementName,
            regionCode as double
        )

        measurementList.close()
    }

    // --------------------------------------------------------
    // MOSTRAR PROGRESO
    // --------------------------------------------------------

    if (
        progressInterval > 0 &&
        (detectionIndex + 1) % progressInterval == 0
    ) {

        double percentage =
            (detectionIndex + 1) *
            100.0 /
            detections.size()

        double elapsedMinutes =
            (
                System.currentTimeMillis() -
                startTime
            ) /
            60000.0

        println(
            String.format(
                'Progreso: %d/%d, %.1f%%, tiempo %.1f min',
                detectionIndex + 1,
                detections.size(),
                percentage,
                elapsedMinutes
            )
        )
    }
}

// ============================================================
// ACTUALIZAR LAS CLASIFICACIONES EN QUPATH
// ============================================================

/*
 * No se llama a resolveHierarchy().
 *
 * Se notifican directamente las clasificaciones modificadas.
 */
if (!changedObjects.isEmpty()) {

    hierarchy.fireObjectClassificationsChangedEvent(
        this,
        changedObjects
    )
}

println ''
println 'Clasificaciones actualizadas en QuPath.'

// ============================================================
// EXPORTAR MEDICIONES A UN ARCHIVO TEMPORAL
// ============================================================

println ''
println '=========================================='
println 'EXPORTACIÓN DE MEDICIONES'
println '=========================================='
println 'Archivo temporal: ' + temporaryMeasurementsPath
println 'Archivo final: ' + detectionMeasurementsPath
println 'Exportando detecciones...'

try {

    /*
     * QuPath genera primero el archivo con los encabezados
     * originales AF405, AF488, AF555, AF647 y AF750.
     */
    saveDetectionMeasurements(
        temporaryMeasurementsPath
    )

} catch (Exception exception) {

    println ''
    println 'ERROR: No se pudieron exportar las mediciones.'
    println 'Mensaje: ' + exception.getMessage()

    throw exception
}

// ============================================================
// VERIFICAR EL ARCHIVO TEMPORAL
// ============================================================

def temporaryFile = new File(
    temporaryMeasurementsPath
)

if (!temporaryFile.exists()) {

    throw new RuntimeException(
        'La exportación terminó, pero no se encontró el ' +
        'archivo temporal: ' +
        temporaryMeasurementsPath
    )
}

if (temporaryFile.length() == 0) {

    throw new RuntimeException(
        'El archivo temporal se creó, pero está vacío: ' +
        temporaryMeasurementsPath
    )
}

// ============================================================
// RENOMBRAR LOS CANALES EN EL ENCABEZADO
// ============================================================

println ''
println 'Renombrando canales en el encabezado...'

def finalOutputFile = new File(
    detectionMeasurementsPath
)

/*
 * Si ya existe un archivo final de una ejecución anterior,
 * será reemplazado.
 */
if (finalOutputFile.exists()) {

    boolean deletedPrevious =
        finalOutputFile.delete()

    if (!deletedPrevious) {

        throw new RuntimeException(
            'No se pudo reemplazar el archivo existente: ' +
            detectionMeasurementsPath
        )
    }
}

int renamedColumns = 0
int totalColumns = 0
long copiedDataRows = 0

/*
 * El archivo se procesa línea por línea.
 * Esto evita cargar todas las detecciones en memoria.
 */
temporaryFile.withReader('UTF-8') { reader ->

    finalOutputFile.withWriter('UTF-8') { writer ->

        String headerLine = reader.readLine()

        if (headerLine == null) {

            throw new RuntimeException(
                'El archivo temporal no contiene encabezado.'
            )
        }

        /*
         * saveDetectionMeasurements genera normalmente
         * un archivo separado por tabuladores.
         */
        String delimiter = '\t'

        def originalHeaders = headerLine.split(
            delimiter,
            -1
        )

        totalColumns = originalHeaders.length

        def renamedHeaders =
            originalHeaders.collect { originalHeader ->
                String renamedHeader = originalHeader

                channelNameMapping.each {
                    oldChannelName,
                    newChannelName ->

                    /*
                     * Formato real exportado por QuPath:
                     *
                     * Nucleus: AF405: Mean
                     * Cytoplasm: AF488: Median
                     * Membrane: AF555: Max
                     * Cell: AF647: Std.Dev.
                     *
                     * El canal está entre separadores ": " y no al
                     * principio del encabezado. Se reemplaza únicamente
                     * un componente completo del encabezado, evitando
                     * coincidencias parciales dentro de otros textos.
                     */
                    String channelAtStart =
                        oldChannelName + ':'

                    String channelInsideHeader =
                        ': ' + oldChannelName + ':'

                    if (renamedHeader.startsWith(channelAtStart)) {
                        renamedHeader =
                            newChannelName +
                            renamedHeader.substring(
                                oldChannelName.length()
                            )
                        renamedColumns++
                    } else if (
                        renamedHeader.contains(channelInsideHeader)
                    ) {
                        renamedHeader = renamedHeader.replace(
                            channelInsideHeader,
                            ': ' + newChannelName + ':'
                        )
                        renamedColumns++
                    }
                }

                return renamedHeader
            }
        writer.write(
            renamedHeaders.join(delimiter)
        )

        writer.newLine()

        /*
         * Copiar todas las filas de datos sin alterar
         * los valores numéricos.
         */
        String dataLine

        while (
            (dataLine = reader.readLine()) != null
        ) {

            writer.write(dataLine)
            writer.newLine()

            copiedDataRows++
        }
    }
}

// ============================================================
// VERIFICAR EL CSV FINAL
// ============================================================

if (!finalOutputFile.exists()) {

    throw new RuntimeException(
        'No se encontró el CSV final: ' +
        detectionMeasurementsPath
    )
}

if (finalOutputFile.length() == 0) {

    throw new RuntimeException(
        'El CSV final está vacío: ' +
        detectionMeasurementsPath
    )
}

// ============================================================
// COMPROBAR LOS ENCABEZADOS FINALES
// ============================================================

String finalHeader = null

finalOutputFile.withReader('UTF-8') { reader ->
    finalHeader = reader.readLine()
}

if (finalHeader == null) {

    throw new RuntimeException(
        'No se pudo leer el encabezado del CSV final.'
    )
}

def oldChannelsStillPresent = []
def newChannelsNotFound = []

channelNameMapping.each {
    oldChannelName,
    newChannelName ->

    boolean oldChannelPresent =
        finalHeader.startsWith(oldChannelName + ':') ||
        finalHeader.contains('\t' + oldChannelName + ':') ||
        finalHeader.contains(': ' + oldChannelName + ':')

    boolean newChannelPresent =
        finalHeader.startsWith(newChannelName + ':') ||
        finalHeader.contains('\t' + newChannelName + ':') ||
        finalHeader.contains(': ' + newChannelName + ':')

    if (oldChannelPresent) {
        oldChannelsStillPresent.add(oldChannelName)
    }

    if (!newChannelPresent) {
        newChannelsNotFound.add(newChannelName)
    }
}

if (!oldChannelsStillPresent.isEmpty()) {

    println ''
    println(
        'ADVERTENCIA: Todavía aparecen estos nombres antiguos: ' +
        oldChannelsStillPresent.join(', ')
    )
}

if (!newChannelsNotFound.isEmpty()) {

    println ''
    println(
        'ADVERTENCIA: No se encontraron columnas para: ' +
        newChannelsNotFound.join(', ')
    )
}

// ============================================================
// ELIMINAR O CONSERVAR EL ARCHIVO TEMPORAL
// ============================================================

if (keepOriginalCSV) {

    println ''
    println 'Se conserva el CSV con los nombres originales:'
    println temporaryMeasurementsPath

} else {

    boolean temporaryDeleted =
        temporaryFile.delete()

    if (temporaryDeleted) {

        println ''
        println 'Archivo temporal eliminado correctamente.'

    } else {

        println ''
        println 'ADVERTENCIA: No se pudo eliminar el archivo:'
        println temporaryMeasurementsPath
    }
}

// ============================================================
// INFORME FINAL
// ============================================================

double totalMinutes =
    (
        System.currentTimeMillis() -
        startTime
    ) /
    60000.0

double fileSizeMegabytes =
    finalOutputFile.length() /
    1024.0 /
    1024.0

println ''
println '=========================================='
println 'CLASIFICACIÓN Y EXPORTACIÓN COMPLETADAS'
println '=========================================='
println 'Imagen: ' + imageName
println 'Detecciones totales: ' + detections.size()

println(
    'Filas exportadas: ' +
    copiedDataRows
)

println(
    'Clasificaciones modificadas: ' +
    changedClassifications
)

println(
    'Clasificaciones sin cambios: ' +
    unchangedClassifications
)

println(
    'Detecciones fuera de anotaciones: ' +
    outsideAnnotations
)

println(
    'Detecciones con ROI no válido: ' +
    invalidDetections
)

println(
    'Errores de contención: ' +
    containmentErrors
)

println ''
println 'Recuento por Classification:'

annotationClassPriority.eachWithIndex {
    classDefinition,
    index ->

    String detectionClassName =
        classDefinition.detectionClassName as String

    println(
        '  ' +
        detectionClassName +
        ': ' +
        classificationCounts[detectionClassName] +
        ' | código ' +
        (index + 1)
    )
}

println ''
println 'Canales renombrados en el CSV:'

channelNameMapping.each {
    oldChannelName,
    newChannelName ->

    println(
        '  ' +
        oldChannelName +
        ' -> ' +
        newChannelName
    )
}

println ''
println 'Columnas totales: ' + totalColumns

println(
    'Columnas de intensidad renombradas: ' +
    renamedColumns
)

println 'CSV final: ' + detectionMeasurementsPath

println(
    String.format(
        'Tamaño del CSV: %.2f MB',
        fileSizeMegabytes
    )
)

println(
    String.format(
        'Tiempo total: %.2f minutos',
        totalMinutes
    )
)

println ''
println 'No se ejecutó resolveHierarchy().'
println 'La columna Classification se asignó directamente.'

println ''
println 'Valores de Classification:'
println '  Annotation (Lesion_dermis)'
println '  Annotation (Lesion_epidermis)'
println '  Annotation (Lesion)'
println '  Annotation (dermis)'
println '  Annotation (epidermis)'

println ''
println 'Nombres finales de canales:'
println '  DAPI'
println '  K14'
println '  MDK'
println '  PTN'
println '  S100B'
println '=========================================='
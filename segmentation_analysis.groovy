import qupath.lib.objects.PathObject
import qupath.lib.objects.PathObjects
import qupath.lib.objects.classes.PathClassFactory
import qupath.lib.images.ImageData
import qupath.lib.gui.scripting.QPEx

// Import Cellpose extension
import qupath.ext.biop.cellpose.Cellpose2D

// =======================================================
// CONFIGURACIÓN
// =======================================================

// Nombre del pixel classifier guardado en QuPath
def pixelClassifierName = "model_442"

// Modelo de Cellpose
def cellposeModel = "model_nuclei"

// Carpeta donde se exportarán las measurements
def outputDir = buildFilePath(PROJECT_BASE_DIR, "measurements_cellpose")
mkdirs(outputDir)

// Parámetros para crear anotaciones desde pixel classifier
def minArea = 10000
def minHoleArea = 1000

// Clases que queremos eliminar después del pixel classification
def classesToRemove = ["BG", "foliculo", "follicle", "Follicle"]

// Clases que queremos conservar y sobre las que correr Cellpose
def classesForCellpose = ["dermis", "epidermis"]

// Primer plano Z.
// En QuPath normalmente el primer Z es 0.
def firstZ = 0

// Parámetros Cellpose
def cellposeChannel = "AF405"
def cellposePixelSize = 0.340
def cellposeDiameter = 20.0
def cellprobThreshold = 1.8
def flowThreshold = 0.8
def cellExpansion = 3.0

// =======================================================
// OBTENER IMAGEN ACTUAL
// =======================================================

def imageData = getCurrentImageData()

if (imageData == null) {
    println "ERROR: no hay ninguna imagen abierta."
    return
}

def hierarchy = imageData.getHierarchy()
def entry = getProjectEntry()

def imageName

if (entry != null) {
    imageName = entry.getImageName()
} else {
    imageName = imageData.getServer().getMetadata().getName()
}

println "======================================================="
println "Procesando únicamente la imagen abierta:"
println imageName
println "======================================================="



// =======================================================
// 2. EJECUTAR PIXEL CLASSIFIER Y CREAR ANOTACIONES
// ======================================================


def annotationsAfterClassifier = getAnnotationObjects()

println "Anotaciones creadas por pixel classifier: " + annotationsAfterClassifier.size()

// =======================================================
// 3. BORRAR BG, FOLICULO Y CLASES NO DESEADAS
// =======================================================

println "Filtrando clases: eliminando " + classesToRemove

def annotationsToDelete = []

annotationsAfterClassifier.each { ann ->

    def pathClass = ann.getPathClass()
    def className = pathClass == null ? "" : pathClass.toString()

    if (classesToRemove.contains(className)) {
        annotationsToDelete << ann
    }
}

if (!annotationsToDelete.isEmpty()) {
    println "Eliminando anotaciones BG/foliculo: " + annotationsToDelete.size()
    hierarchy.removeObjects(annotationsToDelete, true)
} else {
    println "No se encontraron anotaciones BG/foliculo para eliminar."
}

fireHierarchyUpdate()

// =======================================================
// 4. MANTENER SOLO DERMIS Y EPIDERMIS
// =======================================================

println "Manteniendo solo clases: " + classesForCellpose

def annotationsCurrent = getAnnotationObjects()
def annotationsNotWanted = []

annotationsCurrent.each { ann ->

    def pathClass = ann.getPathClass()
    def className = pathClass == null ? "" : pathClass.toString()

    if (!classesForCellpose.contains(className)) {
        annotationsNotWanted << ann
    }
}

if (!annotationsNotWanted.isEmpty()) {
    println "Eliminando anotaciones de otras clases: " + annotationsNotWanted.size()
    hierarchy.removeObjects(annotationsNotWanted, true)
} else {
    println "Todas las anotaciones restantes son dermis/epidermis."
}

fireHierarchyUpdate()

// =======================================================
// 5. MANTENER SOLO ANOTACIONES DEL PRIMER Z
// =======================================================

println "Filtrando anotaciones para conservar solo Z = " + firstZ

def annotationsForZFiltering = getAnnotationObjects()
def annotationsOtherZ = []

annotationsForZFiltering.each { ann ->

    def roi = ann.getROI()

    if (roi == null) {
        annotationsOtherZ << ann
        return
    }

    def z = roi.getZ()

    if (z != firstZ) {
        annotationsOtherZ << ann
    }
}

if (!annotationsOtherZ.isEmpty()) {
    println "Eliminando anotaciones de otros planos Z: " + annotationsOtherZ.size()
    hierarchy.removeObjects(annotationsOtherZ, true)
} else {
    println "No se encontraron anotaciones fuera de Z = " + firstZ
}

fireHierarchyUpdate()

// =======================================================
// 6. SELECCIONAR DERMIS Y EPIDERMIS FINALES
// =======================================================

def annotationsForCellpose = getAnnotationObjects().findAll { ann ->

    def pathClass = ann.getPathClass()
    def className = pathClass == null ? "" : pathClass.toString()

    return classesForCellpose.contains(className)
}

println "Anotaciones finales para Cellpose: " + annotationsForCellpose.size()

annotationsForCellpose.eachWithIndex { ann, i ->
    def className = ann.getPathClass() == null ? "No class" : ann.getPathClass().toString()
    def z = ann.getROI() == null ? "NA" : ann.getROI().getZ()
    println "  Annotation " + i + " | class = " + className + " | z = " + z
}

// =======================================================
// 7. EJECUTAR CELLPOSE SOBRE DERMIS Y EPIDERMIS
// =======================================================

if (annotationsForCellpose.isEmpty()) {

    println "WARNING: no hay anotaciones dermis/epidermis en Z = " + firstZ
    println "Se salta Cellpose."

} else {

    println "Ejecutando Cellpose..."
    println "Modelo: " + cellposeModel
    println "Canal: " + cellposeChannel
    println "Pixel size: " + cellposePixelSize
    println "Diameter: " + cellposeDiameter
    println "cellprobThreshold: " + cellprobThreshold
    println "flowThreshold: " + flowThreshold
    println "cellExpansion: " + cellExpansion

    selectObjects(annotationsForCellpose)

    def cellpose = Cellpose2D.builder(cellposeModel)
            .pixelSize(cellposePixelSize)
            .channels(cellposeChannel)
            .diameter(cellposeDiameter)
            //.measureShape()
            .measureIntensity()
            .cellprobThreshold(cellprobThreshold)
            .flowThreshold(flowThreshold)
            .cellExpansion(cellExpansion)
            .build()

    cellpose.detectObjects(imageData, annotationsForCellpose)

    fireHierarchyUpdate()

    println "Cellpose terminado."
}

// =======================================================
// 8. EXPORTAR MEASUREMENTS
// =======================================================

println "Exportando measurements..."

def safeName = imageName
        .replaceAll("[^a-zA-Z0-9_\\-\\.]", "_")
        .replaceAll("\\.czi", "")
        .replaceAll("\\.CZI", "")
        .replaceAll("\\.ome\\.tiff", "")
        .replaceAll("\\.ome\\.tif", "")
        .replaceAll("\\.tiff", "")
        .replaceAll("\\.tif", "")

def annotationMeasurementsPath = buildFilePath(
        outputDir,
        safeName + "_annotation.csv"
)

def detectionMeasurementsPath = buildFilePath(
        outputDir,
        safeName + "_single_cell.csv"
)

//saveAnnotationMeasurements(annotationMeasurementsPath)
saveDetectionMeasurements(detectionMeasurementsPath)

println "Measurements de anotaciones guardadas en:"
println annotationMeasurementsPath

println "Measurements de detections/células guardadas en:"
println detectionMeasurementsPath

// =======================================================
// 9. GUARDAR CAMBIOS EN EL PROYECTO
// =======================================================

if (entry != null) {
    entry.saveImageData(imageData)
    println "Cambios guardados en el proyecto."
} else {
    println "WARNING: la imagen no está asociada a una entrada de proyecto."
    println "Los resultados existen en la sesión actual, pero puede que no se guarden en el proyecto."
}

println "======================================================="
println "Procesamiento terminado para imagen:"
println imageName
println "Resultados exportados en:"
println outputDir
println "======================================================="

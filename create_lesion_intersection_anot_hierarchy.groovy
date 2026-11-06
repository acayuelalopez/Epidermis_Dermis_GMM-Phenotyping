import qupath.lib.objects.PathObjects
import qupath.lib.roi.GeometryTools

// ============================================================
// CONFIGURACIÓN
// ============================================================

// Clases originales
String dermisClassName = 'dermis'
String epidermisClassName = 'epidermis'
String lesionClassName = 'Lesion'

// Clases nuevas para las intersecciones
String lesionDermisClassName = 'Lesion_dermis'
String lesionEpidermisClassName = 'Lesion_epidermis'

// Eliminar intersecciones anteriores antes de recalcularlas.
// Esto evita duplicados al ejecutar el script varias veces.
boolean removePreviousIntersections = true

// Resolver la jerarquía al terminar.
boolean resolveObjectHierarchy = true

// ============================================================
// OBTENER IMAGEN Y JERARQUÍA
// ============================================================

def imageData = getCurrentImageData()

if (imageData == null) {
    println 'ERROR: No hay ninguna imagen abierta.'
    return
}

def hierarchy = imageData.getHierarchy()
def annotations = getAnnotationObjects()

if (annotations == null || annotations.isEmpty()) {
    println 'ERROR: La imagen no contiene anotaciones.'
    return
}

// ============================================================
// FUNCIÓN PARA COMPARAR CLASES
// ============================================================

/*
 * En la interfaz de QuPath se muestra:
 *
 * Annotation (dermis)
 * Annotation (epidermis)
 * Annotation (Lesion)
 *
 * Sin embargo, el nombre real de la clase es:
 *
 * dermis
 * epidermis
 * Lesion
 */
def hasClassName = { annotation, String expectedName ->

    def pathClass = annotation.getPathClass()

    if (pathClass == null) {
        return false
    }

    return pathClass.toString().equalsIgnoreCase(expectedName)
}

// ============================================================
// MOSTRAR INFORMACIÓN INICIAL
// ============================================================

println ''
println '=========================================='
println 'GENERACIÓN DE INTERSECCIONES'
println '=========================================='
println 'Imagen: ' + imageData.getServer().getMetadata().getName()
println 'Anotaciones iniciales: ' + annotations.size()
println ''

// ============================================================
// ELIMINAR INTERSECCIONES CREADAS ANTERIORMENTE
// ============================================================

int removedPrevious = 0

if (removePreviousIntersections) {

    def previousIntersections = annotations.findAll { annotation ->

        hasClassName(annotation, lesionDermisClassName) ||
        hasClassName(annotation, lesionEpidermisClassName)
    }

    if (!previousIntersections.isEmpty()) {

        println(
            'Eliminando intersecciones anteriores: ' +
            previousIntersections.size()
        )

        /*
         * false evita eliminar los objetos hijos.
         * Después se resolverá nuevamente toda la jerarquía.
         */
        hierarchy.removeObjects(
            previousIntersections,
            false
        )

        removedPrevious = previousIntersections.size()

    } else {

        println 'No se encontraron intersecciones anteriores.'
    }

    // Actualizar la lista después de eliminar objetos
    annotations = getAnnotationObjects()
}

// ============================================================
// LOCALIZAR LAS ANOTACIONES ORIGINALES
// ============================================================

def dermisAnnotations = annotations.findAll { annotation ->
    hasClassName(annotation, dermisClassName)
}

def epidermisAnnotations = annotations.findAll { annotation ->
    hasClassName(annotation, epidermisClassName)
}

def lesionAnnotations = annotations.findAll { annotation ->
    hasClassName(annotation, lesionClassName)
}

// ============================================================
// MOSTRAR ANOTACIONES LOCALIZADAS
// ============================================================

println ''
println 'Anotaciones encontradas:'
println '  Annotation (dermis): ' + dermisAnnotations.size()
println '  Annotation (epidermis): ' + epidermisAnnotations.size()
println '  Annotation (Lesion): ' + lesionAnnotations.size()
println ''

// ============================================================
// FUNCIÓN PARA RESOLVER LA JERARQUÍA
// ============================================================

def resolveHierarchySafely = {

    if (!resolveObjectHierarchy) {
        println 'Resolución de jerarquía desactivada.'
        return
    }

    println ''
    println 'Resolviendo la jerarquía de objetos...'
    println 'Este proceso puede tardar si hay muchas detecciones.'

    try {

        hierarchy.resolveHierarchy()

        println 'Jerarquía resuelta correctamente.'

    } catch (Exception exception) {

        println ''
        println 'ADVERTENCIA: No se pudo resolver la jerarquía.'
        println 'Mensaje: ' + exception.getMessage()
    }

    hierarchy.fireHierarchyChangedEvent(this)
}

// ============================================================
// SI NO EXISTE LESION, NO GENERAR NADA
// ============================================================

if (lesionAnnotations.isEmpty()) {

    println 'No existe ninguna Annotation (Lesion).'
    println 'No se crearán anotaciones de intersección.'
    println ''
    println 'Se conservarán las anotaciones existentes:'
    println '  Annotation (dermis): ' + dermisAnnotations.size()
    println '  Annotation (epidermis): ' + epidermisAnnotations.size()

    // Resolver la jerarquía aunque no exista lesión
    resolveHierarchySafely()

    println ''
    println '=========================================='
    println 'PROCESAMIENTO FINALIZADO'
    println '=========================================='

    return
}

// ============================================================
// COMPROBAR SI EXISTEN DERMIS O EPIDERMIS
// ============================================================

if (dermisAnnotations.isEmpty()) {
    println 'ADVERTENCIA: No existe Annotation (dermis).'
    println 'No se podrá crear Annotation (Lesion_dermis).'
}

if (epidermisAnnotations.isEmpty()) {
    println 'ADVERTENCIA: No existe Annotation (epidermis).'
    println 'No se podrá crear Annotation (Lesion_epidermis).'
}

// ============================================================
// OBTENER O CREAR LAS CLASES DE INTERSECCIÓN
// ============================================================

def lesionDermisClass = getPathClass(
    lesionDermisClassName
)

def lesionEpidermisClass = getPathClass(
    lesionEpidermisClassName
)

// ============================================================
// PREPARAR CONTADORES
// ============================================================

def newAnnotations = []

int lesionDermisCreated = 0
int lesionEpidermisCreated = 0

int emptyDermisIntersections = 0
int emptyEpidermisIntersections = 0

int differentPlaneSkipped = 0
int geometryErrors = 0

// ============================================================
// CREAR INTERSECCIONES LESION-DERMIS
// ============================================================

println ''
println 'Calculando intersecciones Lesion-dermis...'

lesionAnnotations.each { lesionAnnotation ->

    def lesionROI = lesionAnnotation.getROI()

    if (lesionROI == null || lesionROI.isEmpty()) {

        println 'ADVERTENCIA: Se encontró una lesión con ROI vacío.'
        return
    }

    dermisAnnotations.each { dermisAnnotation ->

        def dermisROI = dermisAnnotation.getROI()

        if (dermisROI == null || dermisROI.isEmpty()) {
            return
        }

        /*
         * Solo calcular la intersección cuando ambas anotaciones
         * pertenecen al mismo plano Z y al mismo tiempo T.
         */
        if (
            lesionROI.getZ() != dermisROI.getZ() ||
            lesionROI.getT() != dermisROI.getT()
        ) {

            differentPlaneSkipped++
            return
        }

        /*
         * Comprobación rápida de los límites.
         * Evita calcular geometrías si los bounding boxes
         * no se solapan.
         */
        boolean boundingBoxesOverlap =
            lesionROI.getBoundsX() <
                dermisROI.getBoundsX() + dermisROI.getBoundsWidth() &&
            lesionROI.getBoundsX() + lesionROI.getBoundsWidth() >
                dermisROI.getBoundsX() &&
            lesionROI.getBoundsY() <
                dermisROI.getBoundsY() + dermisROI.getBoundsHeight() &&
            lesionROI.getBoundsY() + lesionROI.getBoundsHeight() >
                dermisROI.getBoundsY()

        if (!boundingBoxesOverlap) {
            emptyDermisIntersections++
            return
        }

        try {

            def lesionGeometry = lesionROI.getGeometry()
            def dermisGeometry = dermisROI.getGeometry()

            /*
             * Calcular:
             *
             * Lesion ∩ dermis
             */
            def intersectionGeometry =
                lesionGeometry.intersection(dermisGeometry)

            /*
             * Conservar únicamente componentes con área.
             * Se eliminan puntos y líneas.
             */
            intersectionGeometry =
                GeometryTools.ensurePolygonal(
                    intersectionGeometry
                )

            if (
                intersectionGeometry == null ||
                intersectionGeometry.isEmpty() ||
                intersectionGeometry.getArea() <= 0
            ) {

                emptyDermisIntersections++
                return
            }

            /*
             * Convertir la geometría resultante nuevamente
             * a un ROI de QuPath.
             */
            def intersectionROI =
                GeometryTools.geometryToROI(
                    intersectionGeometry,
                    lesionROI.getImagePlane()
                )

            if (
                intersectionROI == null ||
                intersectionROI.isEmpty()
            ) {

                emptyDermisIntersections++
                return
            }

            /*
             * Crear:
             *
             * Annotation (Lesion_dermis)
             */
            def intersectionAnnotation =
                PathObjects.createAnnotationObject(
                    intersectionROI,
                    lesionDermisClass
                )

            newAnnotations.add(intersectionAnnotation)
            lesionDermisCreated++

        } catch (Exception exception) {

            geometryErrors++

            println(
                'ADVERTENCIA: No se pudo calcular una ' +
                'intersección Lesion-dermis.'
            )

            println 'Mensaje: ' + exception.getMessage()
        }
    }
}

// ============================================================
// CREAR INTERSECCIONES LESION-EPIDERMIS
// ============================================================

println 'Calculando intersecciones Lesion-epidermis...'

lesionAnnotations.each { lesionAnnotation ->

    def lesionROI = lesionAnnotation.getROI()

    if (lesionROI == null || lesionROI.isEmpty()) {
        return
    }

    epidermisAnnotations.each { epidermisAnnotation ->

        def epidermisROI = epidermisAnnotation.getROI()

        if (
            epidermisROI == null ||
            epidermisROI.isEmpty()
        ) {
            return
        }

        /*
         * Solo calcular la intersección cuando las anotaciones
         * pertenecen al mismo plano Z y tiempo T.
         */
        if (
            lesionROI.getZ() != epidermisROI.getZ() ||
            lesionROI.getT() != epidermisROI.getT()
        ) {

            differentPlaneSkipped++
            return
        }

        /*
         * Comprobación rápida de solapamiento entre
         * bounding boxes.
         */
        boolean boundingBoxesOverlap =
            lesionROI.getBoundsX() <
                epidermisROI.getBoundsX() +
                epidermisROI.getBoundsWidth() &&
            lesionROI.getBoundsX() +
                lesionROI.getBoundsWidth() >
                epidermisROI.getBoundsX() &&
            lesionROI.getBoundsY() <
                epidermisROI.getBoundsY() +
                epidermisROI.getBoundsHeight() &&
            lesionROI.getBoundsY() +
                lesionROI.getBoundsHeight() >
                epidermisROI.getBoundsY()

        if (!boundingBoxesOverlap) {
            emptyEpidermisIntersections++
            return
        }

        try {

            def lesionGeometry = lesionROI.getGeometry()
            def epidermisGeometry =
                epidermisROI.getGeometry()

            /*
             * Calcular:
             *
             * Lesion ∩ epidermis
             */
            def intersectionGeometry =
                lesionGeometry.intersection(
                    epidermisGeometry
                )

            /*
             * Eliminar componentes que no tengan área,
             * como líneas o puntos.
             */
            intersectionGeometry =
                GeometryTools.ensurePolygonal(
                    intersectionGeometry
                )

            if (
                intersectionGeometry == null ||
                intersectionGeometry.isEmpty() ||
                intersectionGeometry.getArea() <= 0
            ) {

                emptyEpidermisIntersections++
                return
            }

            /*
             * Convertir el resultado nuevamente a ROI.
             */
            def intersectionROI =
                GeometryTools.geometryToROI(
                    intersectionGeometry,
                    lesionROI.getImagePlane()
                )

            if (
                intersectionROI == null ||
                intersectionROI.isEmpty()
            ) {

                emptyEpidermisIntersections++
                return
            }

            /*
             * Crear:
             *
             * Annotation (Lesion_epidermis)
             */
            def intersectionAnnotation =
                PathObjects.createAnnotationObject(
                    intersectionROI,
                    lesionEpidermisClass
                )

            newAnnotations.add(intersectionAnnotation)
            lesionEpidermisCreated++

        } catch (Exception exception) {

            geometryErrors++

            println(
                'ADVERTENCIA: No se pudo calcular una ' +
                'intersección Lesion-epidermis.'
            )

            println 'Mensaje: ' + exception.getMessage()
        }
    }
}

// ============================================================
// AÑADIR LAS NUEVAS ANOTACIONES
// ============================================================

println ''

if (!newAnnotations.isEmpty()) {

    println(
        'Añadiendo nuevas anotaciones: ' +
        newAnnotations.size()
    )

    hierarchy.addObjects(newAnnotations)

    println 'Nuevas anotaciones añadidas correctamente.'

} else {

    println 'No se generaron nuevas anotaciones.'
    println(
        'La lesión no se solapa con dermis o epidermis, ' +
        'o las anotaciones están en planos diferentes.'
    )
}

// ============================================================
// RESOLVER LA JERARQUÍA
// ============================================================

resolveHierarchySafely()

// ============================================================
// MOSTRAR RESULTADO FINAL
// ============================================================

println ''
println '=========================================='
println 'INTERSECCIONES COMPLETADAS'
println '=========================================='

println(
    'Annotation (dermis) conservadas: ' +
    dermisAnnotations.size()
)

println(
    'Annotation (epidermis) conservadas: ' +
    epidermisAnnotations.size()
)

println(
    'Annotation (Lesion) conservadas: ' +
    lesionAnnotations.size()
)

println(
    'Annotation (Lesion_dermis) creadas: ' +
    lesionDermisCreated
)

println(
    'Annotation (Lesion_epidermis) creadas: ' +
    lesionEpidermisCreated
)

println(
    'Intersecciones Lesion-dermis vacías: ' +
    emptyDermisIntersections
)

println(
    'Intersecciones Lesion-epidermis vacías: ' +
    emptyEpidermisIntersections
)

println(
    'Combinaciones omitidas por plano Z/T diferente: ' +
    differentPlaneSkipped
)

println(
    'Errores geométricos: ' +
    geometryErrors
)

println(
    'Intersecciones anteriores eliminadas: ' +
    removedPrevious
)

println(
    'Total de nuevas anotaciones: ' +
    newAnnotations.size()
)

println 'Jerarquía resuelta: ' + resolveObjectHierarchy
println '=========================================='

println ''
println 'Resultado esperado:'
println '  Annotation (dermis)'
println '  Annotation (epidermis)'
println '  Annotation (Lesion)'

if (lesionDermisCreated > 0) {
    println '  Annotation (Lesion_dermis)'
}

if (lesionEpidermisCreated > 0) {
    println '  Annotation (Lesion_epidermis)'
}

println ''
println 'Recuerda guardar los cambios del proyecto.'
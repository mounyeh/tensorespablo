#===============================================================================
#
#===============================================================================
# Fichero:        MarcacionElipsoide.tcl
# Procesos:  
#   MarcacionElipsoideInit
#   MarcacionElipsoideBuildGUI
#   MarcacionElipsoideBuildHelpFrame
#   MarcacionElipsoideBuildVolumenesFrame
#   MarcacionElipsoideBuildMarcarFrame
#   MarcacionElipsoideBuildParametrosFrame
#   MarcacionElipsoideBuildTrazarFrame
#   MarcacionElipsoideEnter
#   MarcacionElipsoideExit
#   MarcacionElipsoideUpdateGUI
#   MarcacionElipsoideSetOriginal
#   MarcacionElipsoideSetWorking
#   MarcacionElipsoideGetOriginalID
#   MarcacionElipsoideGetWorkingID
#   MarcacionElipsoideCount
#   MarcacionElipsoideShowFile
#   MarcacionElipsoideBindingCallback
#   MarcacionElipsoidePrepareResult
#   MarcacionElipsoidePrepareResultVolume
#   RunMarcacionElipsoide
#   MarcacionElipsoideCarga
#   MarcacionElipsoideTrazar
#   MarcacionElipsoideCalcular
#   MarcacionElipsoideModeloVerosimilitudCalc
#   MarcacionElipsoideModeloOptimizacionSA
#   MarcacionElipsoideEstablece
#   MarcacionElipsoideSetupBeforeApplyEffect
#   MarcacionElipsoideUpdateAfterApplyEffect
#==========================================================================auto=

#-------------------------------------------------------------------------------
#  Descripción
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#  Variables
#  Estas son (algunas de) las variables definidas por este módulo.
# 
#
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# .PROC MarcacionElipsoideInit
#  El proceso "Init" es llamado automáticamente por el slicer.
#  Sitúa la información sobre el módulo en un array global llamado Module, 
#  y también inicializa variables module-level.
# .ARGS
# .END
#-------------------------------------------------------------------------------
proc MarcacionElipsoideInit {} {

    global MarcacionElipsoide Module Volume Model
    
    puts "MarcacionElipsoideInit begin"

    # Definición de Tabs
    #------------------------------------
    # Descripción:
    #   A cada módulo se le da un botón en el menú main del Slicer.
    #   Cuando se presiona ese botón aparece una fila de tabs, y hay un panel
    #   para interfaz con el usuario en cada tab. Si todas las tabs no caben en una
    #   fila, entonces la última tab se crea automáticamente para decir "More", y 
    #   presionarla muestra una segunda fila de tabs.
    #
    #   Define tus tabs aquí como se muestra más abajo. Las opciones son:
    #   
    #   row1List = lista de identificadores las tabs. (los identificadores deben ser una sola palabra que identifique
    #			 unívocamente)
    #   row1Name = lista de nombres para las tabs. (los nombres aparecen en el interfaz con el usuario
    #              y pueden ser no unívocos, con múltiples palabras)
    #   row1,tab = identificador de la tab inicial
    #   row2List = una segunda fila opcional de tabs si la primera fila es demasiado pequeña
    #   row2Name = como en row1
    #   row2,tab = como en row1 
    #
    set m MarcacionElipsoide
    set Module($m,row1List) "Help Volumenes Marcar Parametros Calcular"
    set Module($m,row1Name) "{Help} {Volumenes} {Marcar} {Parametros} {Calcular}"
    set Module($m,row1,tab) Volumenes

    # Definir Procedimientos
    #------------------------------------
    # Descripción:
    #   El Slicer toma como fuentes todos los ficheros *.tcl, y entonces llama a las funciones
    #   Init de cada módulo, seguido por las funciones VTK, y finalmente
    #   las funciones GUI. Una función MRML es llamada allá donde existan cambios en el
    #   árbol MRML debidos a la creación/borrado de nodos.
    #   
    #   Así como el procedimiento Init se requiere para cada módulo, los otros 
    #   procedimientos son opcionales. Si existen, entonces su nombre (el cual
    #   puede ser cualquier cosa) se registra con una línea como esta:
    #
    #   set Module($m,procVTK) MarcacionElipsoideBuildVTK
    #
    #   Las opciones son:
    #
    #   procGUI   = Construye la interfaz gráfica de usuario
    #   procVTK   = Construye los objetos VTK
    #   procMRML  = Actualiza después de cambios en el árbol MRML debidos a la creación
    #               o borrado de nodos.
    #   procEnter = Llamado cuando el usuario entra en este módulo haciendo click en
    #               su botón en el menú principal
    #   procExit  = Se llama cuando el usuario abandona este módulo haciendo click en
    #               el botón de otros módulos
    #   procCameraMotion = Se llama justo antes de que la cámara del renderer 
    #                      activo esté a punto de moverse 
    #   procStorePresets  = Se llama cuando el usuario mantiene pulsado uno de los botones de Presets.
    #   procRecallPresets  = Se llama cuando el usuario pulsa uno de los botones de Presets.
    #               
    #   Note: si usas presets, debes estar seguro que has dado una cadena de preset
    #   por defecto en tu función init, con la forma: 
    #   set Module($m,presets) "key1='val1' key2='val2' ..."
    #   
    set Module($m,procGUI)   MarcacionElipsoideBuildGUI
    set Module($m,procEnter) MarcacionElipsoideEnter
    set Module($m,procExit)  MarcacionElipsoideExit
    set Module($m,procMRML)  MarcacionElipsoideUpdateGUI
    set Module($m,procVTK)   MarcacionElipsoideBuildVTK

    # Definir Dependencias
    #------------------------------------
    # Descripcion:
    #   Grabar cualquier otro módulo del que dependa éste. Se usa 
    #   para comprobar que todos los módulos necesarios se y han cargado cuando el Slicer se ejecuta.
    #
    #Definición de dependencias de nuestro módulo respecto al Fiducials
    set Module($m,depend) Fiducials

    # MarcacionElipsoideEstablecer la información de versión
    #------------------------------------
    # Descripcion:
    #   Grabar el número de versión para grabar ante Help->Version Info.
    #   Las cadenas con el símbolo $ le dicen a CVS que inserte automáticamente el
    #   número de revisión apropiada y la fecha cuando el módulo sea comprobado.
    #   
    lappend Module(versions) [ParseCVSInfo $m \
        {$Revision: 1.0 $} {$Date: 2004/07/26 22:19:00 $}]

    # Inicialización de las variables module-level
    #------------------------------------
    # Descripcion:
    #   Manten un array global con el mismo nombre que el módulo.
    #   Este es un método manejable para organizar las variables globales que
    #   necesitan acceder los procesos en este módulo y en otros.
    #
    set MarcacionElipsoide(count) 0
    set MarcacionElipsoide(InputVol) $Volume(idNone)
    set MarcacionElipsoide(ResultVol) $Model(idNone)
    set MarcacionElipsoide(idOriginal) $Volume(idNone)
    set MarcacionElipsoide(idWorking) NEW
    set MarcacionElipsoide(FileName)  ""
    set MarcacionElipsoide(prefixWorking) ""
    set MarcacionElipsoide(nameWorking) Working
    set MarcacionElipsoide(ResTheta) 32
    set MarcacionElipsoide(ResPhi) 32
    set MarcacionElipsoide(J) 50
    set MarcacionElipsoide(K) 15
    set MarcacionElipsoide(drmax) 20.0
    set MarcacionElipsoide(Ng) 15
    set MarcacionElipsoide(Ns) 25
    set MarcacionElipsoide(Ve0) 35
    set MarcacionElipsoide(Ve1) 35
    set MarcacionElipsoide(Ve2) 35
    set MarcacionElipsoide(Ve3) 35
    set MarcacionElipsoide(Ve4) 35
    set MarcacionElipsoide(Ve5) 35
    set MarcacionElipsoide(Ve6) 35
    set MarcacionElipsoide(AlfaPot0) 30
    set MarcacionElipsoide(AlfaPot1) 30
    set MarcacionElipsoide(AlfaPot2) 10
    set MarcacionElipsoide(AlfaPot3) 10
    set MarcacionElipsoide(AlfaPot4) 10
    set MarcacionElipsoide(AlfaPot5) 4
    set MarcacionElipsoide(AlfaPot6) 6
    set MarcacionElipsoide(BetaPot0) 0.1
    set MarcacionElipsoide(BetaPot1) 0.1
    set MarcacionElipsoide(BetaPot2) 0.1
    set MarcacionElipsoide(BetaPot3) 0.1
    set MarcacionElipsoide(BetaPot4) 0.05
    set MarcacionElipsoide(color) 1
    set MarcacionElipsoide(FactConv) 1
#indica que por defecto se carga la marcación inicial de un fichero "PuntosMarcados.vtk"
    set MarcacionElipsoide(Cargar) 1
    set MarcacionElipsoide(CambiaFOV) 1
    set MarcacionElipsoide(Nombre) "Legibles/RINON1b.VOL"
    
 

    #Lista predefinida para nuestro módulo
    set MarcacionElipsoide(fiducialsListCreated) 0


 
    # ¡Enlaces de eventos! (ver MarcacionElipsoideEnter, MarcacionElipsoideExit, tcl-shared/Events.tcl)
    set MarcacionElipsoide(eventManager)  { \
        {all <Shift-1> {MarcacionElipsoideBindingCallback Shift-1 %W %X %Y %x %y %t}} \
        {all <Shift-2> {MarcacionElipsoideBindingCallback Shift-2 %W %X %Y %x %y %t}} \
        {all <Shift-3> {MarcacionElipsoideBindingCallback Shift-3 %W %X %Y %x %y %t}} }
    
    puts "MarcacionElipsoideInit end"

}


#-------------------------------------------------------------------------------
# .PROC MarcacionElipsoideBuildVTK
# 
# Este proceso construye los objetos VTK necesarios.
# .ARGS
# .END
#-------------------------------------------------------------------------------
proc MarcacionElipsoideBuildVTK {} {

    global MarcacionElipsoide
    

}


#NADA CLARO CÓMO FUNCIONA EL SIGUIENTE PROCESO

#-------------------------------------------------------------------------------
# .PROC MarcacionElipsoideUpdateGUI
# 
# A este proceso se le llama para actualizar los botones
# debidos a cuestiones tales como volúmenes o modelos que están siendo añadidos o quitados.
# (Nota: para hacer esto, este proceso debe ser <modulo>procMRML. Justo ahora,
# estos botones están siendo actualizados automáticamente, ya que han sido añadidos
# a listas actualizadas en VolumesUpdateMRML y ModelsUpdateMRML. Así que este procedimiento
# no se usa corrientemente.)
# .ARGS
# .END
#-------------------------------------------------------------------------------
proc MarcacionElipsoideUpdateGUI {} {
    global MarcacionElipsoide Volume Model
    
    #DevUpdateNodeSelectButton Volume MarcacionElipsoide InputVol InputVol DevSelectNode
    #DevUpdateNodeSelectButton Model  MarcacionElipsoide ResultVol ResultVol DevSelectNode
    
    # See if the volume for each menu actually exists.
    # If not, use the None volume
    #
    set n $Volume(idNone)
    if {[lsearch $Volume(idList) $MarcacionElipsoide(idOriginal)] == -1} {
        MarcacionElipsoideSetOriginal $n
    }
    if {$MarcacionElipsoide(idWorking) != "NEW" && \
        [lsearch $Volume(idList) $MarcacionElipsoide(idWorking)] == -1} {
        MarcacionElipsoideSetWorking NEW
    }

    # Original Volume menu
    #---------------------------------------------------------------------------
    set m $MarcacionElipsoide(mOriginal)
    $m delete 0 end
    foreach v $Volume(idList) {
        $m add command -label [Volume($v,node) GetName] -command \
            "MarcacionElipsoideSetOriginal $v; RenderAll"
    }

    # Working Volume menu
    #---------------------------------------------------------------------------
    set m $MarcacionElipsoide(mWorking)
    $m delete 0 end
    set idWorking ""
    foreach v $Volume(idList) {
        if {$v != $Volume(idNone)} {
            $m add command -label [Volume($v,node) GetName] -command \
                "MarcacionElipsoideSetWorking $v; RenderAll"
        }
    }

    # Always add a NEW option
    $m add command -label NEW -command "MarcacionElipsoideSetWorking NEW; RenderAll"

    # Set the working volume
    MarcacionElipsoideSetWorking $MarcacionElipsoide(idWorking)

    # Working Volume name field  (name for the NEW volume to be created)
    #---------------------------------------------------------------------------
    set v $MarcacionElipsoide(idWorking)
    if {$v != "NEW"} {
        set MarcacionElipsoide(nameWorking) [Volume($v,node) GetName]
    } else {
        set MarcacionElipsoide(nameWorking) Working
    }
}




# Convención de nombres:
#-------------------------------------------------------------------------------
#
# Usar las siguientes letras iniciales para los nombres:
# t  = toplevel
# f  = frame
# mb = menubutton
# m  = menu
# b  = button
# l  = label
# s  = slider
# i  = image
# c  = checkbox
# r  = radiobutton
# e  = entry
#
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# .PROC MarcacionElipsoideBuildGUI
#
# Crea la Interfaz Gráfica de Usuario.
# .END
#-------------------------------------------------------------------------------
proc MarcacionElipsoideBuildGUI {} {
    
    # Un marco se ha construido ya automáticamente para cada tab.
    # Un marco llamado "Parameters" puede ser referenciado como sigue:
    #   
    #     $Module(<Module name>,f<Tab name>)
    #
    # ie: $Module(MarcacionElipsoide,fMain)
    
    # Este es un bloque de comentarios útil que hace la lectura de esto fácil para todo:
    #-------------------------------------------
    # Jerarquía de marcos:
    #-------------------------------------------
    # Help
    # Volumenes
    # Marcar
    #-------------------------------------------
    
    MarcacionElipsoideBuildHelpFrame
       
    MarcacionElipsoideBuildVolumenesFrame
       
    MarcacionElipsoideBuildMarcarFrame
    
    MarcacionElipsoideBuildParametrosFrame
    
    MarcacionElipsoideBuildCalcularFrame

}

#-------------------------------------------------------------------------------
# .PROC MarcacionElipsoideBuildHelpFrame
#
#   Crea el marco Help
# .END
#-------------------------------------------------------------------------------
proc MarcacionElipsoideBuildHelpFrame {} {


    #-------------------------------------------
    # Marco Help
    #-------------------------------------------
    
    # Escribe la "ayuda" en la forma de pseudo-html.  
    # Referirse a la documentación para detalles sobre la sintaxis.
    #
    set help "
    El módulo MarcacionElipsoide contiene un generador de elipsoide que aproxima volúmenes
    con esta forma desarrollado por Lucilio Cordero Grande.
    Se marcan seis puntos que determinan los tres ejes principales (consecutivamente)
    y se traza el elipsoide.
    <P>
    Los parámetros de entrada son:
    <BR>
    <UL>
    <LI><B> Input image:</B>
    "
    regsub -all "\n" $help {} help
    MainHelpApplyTags MarcacionElipsoide $help
    MainHelpBuildGUI  MarcacionElipsoide

}
# end MarcacionElipsoideBuildHelpFrame

#-------------------------------------------------------------------------------
# .PROC MarcacionElipsoideBuildVolumenesFrame
#
#   Crea el marco Volumenes
# .END
#-------------------------------------------------------------------------------
proc MarcacionElipsoideBuildVolumenesFrame {} {

    global Gui MarcacionElipsoide Module Volume

    #-------------------------------------------
    # Marco Volumenes
    #-------------------------------------------
    set fVolumenes $Module(MarcacionElipsoide,fVolumenes)
    set f $fVolumenes

    frame $f.fHelp      -bg $Gui(activeWorkspace)
    frame $f.fOriginal  -bg $Gui(activeWorkspace) -relief groove -bd 3
    frame $f.fWorking   -bg $Gui(activeWorkspace) -relief groove -bd 3
    frame $f.fVolFilt   -bg $Gui(activeWorkspace)
    frame $f.fStart     -bg $Gui(activeWorkspace) 

    pack  $f.fHelp $f.fOriginal $f.fWorking $f.fVolFilt $f.fStart\
        -side top -padx $Gui(pad) -pady $Gui(pad) -fill x
    
    #-------------------------------------------
    # Volumenes->Help
    #-------------------------------------------
    set f $fVolumenes.fHelp
    
    eval {label $f.l -text \
        "First choose the volumes to edit.\n Type a name for any NEW labelmap.\nThen click `Start'."} \
        $Gui(WLA)
    pack $f.l

    #-------------------------------------------
    # Volumenes->Original
    #-------------------------------------------
    set f $fVolumenes.fOriginal
    
    frame $f.fMenu -bg $Gui(activeWorkspace)
    
    pack $f.fMenu -side top -pady $Gui(pad) -fill x

    #-------------------------------------------
    # Volumenes->Original->Menu
    #-------------------------------------------
    set f $fVolumenes.fOriginal.fMenu
    
    # Volume menu
    eval {label $f.lOriginal -text "Original Grayscale:"} $Gui(WTA)
    
    eval {menubutton $f.mbOriginal -text "None" -relief raised -bd 2 -width 18 \
        -menu $f.mbOriginal.m} $Gui(WMBA)
    eval {menu $f.mbOriginal.m} $Gui(WMA)
    TooltipAdd $f.mbOriginal "Choose the input grayscale volume for editing."
    pack $f.lOriginal -padx $Gui(pad) -side left -anchor e
    pack $f.mbOriginal -padx $Gui(pad) -side left -anchor w
    
    # Save widgets for changing
    set MarcacionElipsoide(mbOriginal) $f.mbOriginal
    set MarcacionElipsoide(mOriginal)  $f.mbOriginal.m
    
    #-------------------------------------------
    # Volumenes->Working
    #-------------------------------------------
    set f $fVolumenes.fWorking
    
    frame $f.fMenu -bg $Gui(activeWorkspace)
    frame $f.fName -bg $Gui(activeWorkspace)
    
    pack $f.fMenu -side top -pady $Gui(pad)
    pack $f.fName -side top -pady $Gui(pad) -fill x

    #-------------------------------------------
    # Volumenes->Working->Menu
    #-------------------------------------------
    set f $fVolumenes.fWorking.fMenu
    
    # Volume menu
    eval {label $f.lWorking -text "Working Labelmap:"} $Gui(WTA)
    
    eval {menubutton $f.mbWorking -text "NEW" -relief raised -bd 2 -width 18 \
        -menu $f.mbWorking.m} $Gui(WMBA)
    eval {menu $f.mbWorking.m} $Gui(WMA)
    TooltipAdd $f.mbWorking "Choose a labelmap to edit, or NEW for a new one."
    pack $f.lWorking $f.mbWorking -padx $Gui(pad) -side left
    
    # Save widgets for changing
    set MarcacionElipsoide(mbWorking) $f.mbWorking
    set MarcacionElipsoide(mWorking)  $f.mbWorking.m
    
    #-------------------------------------------
    # Volumenes->Working->Name
    #-------------------------------------------
    set f $fVolumenes.fWorking.fName
    
    eval {label $f.l -text "Descriptive Name:"} $Gui(WLA)
    eval {entry $f.e -textvariable MarcacionElipsoide(nameWorking)} $Gui(WEA)
    bind $f.e <Return> "set MarcacionElipsoide(idWorking) MarcacionElipsoideGetWorkingID"
    TooltipAdd $f.e "Nickname your NEW volume."
    pack $f.l -padx 3 -side left
    pack $f.e -padx 3 -side left -expand 1 -fill x
    
    # Save widget for disabling name field if not NEW volume
    set MarcacionElipsoide(eNameWorking) $f.e


   #-------------------------------------------
    # Volumenes->Start
    #-------------------------------------------
    set f $fVolumenes.fVolFilt
    
    DevAddButton $f.bVolFilt "Añadir volumen filtrado" "MarcacionElipsoideVolumenFiltrado"
    TooltipAdd $f.bVolFilt "Go go go!"
    pack $f.bVolFilt -side top -padx $Gui(pad) -pady $Gui(pad)
    
    
    #-------------------------------------------
    # Volumenes->Start
    #-------------------------------------------
    set f $fVolumenes.fStart
    
    DevAddButton $f.bStart "Start Editing" "Tab MarcacionElipsoide row1 Marcar"
    TooltipAdd $f.bStart "Go go go!"  
    pack $f.bStart -side top -padx $Gui(pad) -pady $Gui(pad)
     
    
}
# end MarcacionElipsoideBuildVolumenesFrame



#-------------------------------------------------------------------------------
# .PROC MarcacionElipsoideBuildMarcarFrame
#
#   Crea el marco Marcar
# .END
#-------------------------------------------------------------------------------
proc MarcacionElipsoideBuildMarcarFrame {} {

    global Gui MarcacionElipsoide Module Volume

    #-------------------------------------------
    # Marco Marcar
    #-------------------------------------------
    set fMarcar $Module(MarcacionElipsoide,fMarcar)
    set f $fMarcar

    #Para tener acceso a la lista de puntos sin ir a "Fiducials"
    FiducialsAddActiveListFrame $f 7 25 "reformat"
    
    entry $f.eCreateList -width 15 -textvariable Fiducials(newListName)
    bind $f.eCreateList <Return> {set MarcacionElipsoide(identif) \
      [FiducialsCreateFiducialsList "default" $Fiducials(newListName)]
	FiducialsSetActiveList $Fiducials(newListName)
	vtkMarcacionElipsoide MarcacionElipsoide(Generador,$MarcacionElipsoide(identif))
	vtkPlantillaAjustada MarcacionElipsoide(Plantilla,$MarcacionElipsoide(identif))
	vtkFuncionVerosimilitud MarcacionElipsoide(Verosimil,$MarcacionElipsoide(identif))
	vtkOptimizaContorno MarcacionElipsoide(Optimiza,$MarcacionElipsoide(identif))
	vtkValidacion MarcacionElipsoide(Valida,$MarcacionElipsoide(identif))
	vtkVolusonReader MarcacionElipsoide(Voluson,$MarcacionElipsoide(identif))
	
	}
    pack $f.eCreateList -side top


    #Añadimos un botón para escoger los parámetros del algoritmo
    DevAddButton $f.b "Parametros" "Tab MarcacionElipsoide row1 Parametros"
    TooltipAdd $f.b "Go go go"
    pack $f.b -side bottom -padx $Gui(pad) -pady $Gui(pad)
    
    #-------------------------------------------
    # Volumenes->V
    #-------------------------------------------

		eval {label $f.lREs -text "Nombre del fichero"} $Gui(WLA)
    pack $f.lREs -padx $Gui(pad) -pady $Gui(pad)
    
    entry $f.eREs -width 25 -textvariable MarcacionElipsoide(Nombre)
    bind $f.eREs <Return> "MarcacionElipsoideEstablece 27"
    TooltipAdd $f.eREs "Introduce el path del fichero"
    pack $f.eREs

}
# end MarcacionElipsoideBuildMarcarFrame

#-------------------------------------------------------------------------------
# .PROC MarcacionElipsoideBuildParametrosFrame
#
#   Crea el marco Parametros
# .END
#-------------------------------------------------------------------------------
proc MarcacionElipsoideBuildParametrosFrame {} {

    global Gui MarcacionElipsoide Module Volume

    #-------------------------------------------
    # Marco Parametros
    #-------------------------------------------
    set fParametros $Module(MarcacionElipsoide,fParametros)
    set f $fParametros

    eval {label $f.lRes -text "Resolución del elipsoide (Longitud,Latitud)"} $Gui(WLA)
    pack $f.lRes -side top
    
    frame $f.fRes -width 10
    entry $f.fRes.eTheta -width 3 -textvariable MarcacionElipsoide(ResTheta)
    bind $f.fRes.eTheta <Return> "MarcacionElipsoideEstablece 0"
    pack $f.fRes.eTheta -side left -in $f.fRes
    entry $f.fRes.ePhi -width 3 -textvariable MarcacionElipsoide(ResPhi)
    bind $f.fRes.ePhi <Return> "MarcacionElipsoideEstablece 1"
    pack $f.fRes.ePhi -side left -in $f.fRes
    pack $f.fRes -side top
    
    eval {label $f.lJ -text "Número de rayos"} $Gui(WLA)
    pack $f.lJ -side top
    
    entry $f.eJ -width 5 -textvariable MarcacionElipsoide(J)
    bind $f.eJ <Return> "MarcacionElipsoideEstablece 2"
    pack $f.eJ -side top

    eval {label $f.lK -text "Número de puntos de deformación por rayo"} $Gui(WLA)
    pack $f.lK -side top

    entry $f.eK -width 5 -textvariable MarcacionElipsoide(K)
    bind $f.eK <Return> "MarcacionElipsoideEstablece 3"
    pack $f.eK -side top
    
    eval {label $f.ldrmax -text "Radio máximo de la deformación (RAS)"} $Gui(WLA)
    pack $f.ldrmax -side top

    entry $f.edrmax -width 5 -textvariable MarcacionElipsoide(drmax)
    bind $f.edrmax <Return> "MarcacionElipsoideEstablece 4"
    pack $f.edrmax -side top

    eval {label $f.lNg -text "Diámetro del suavizado gaussiano"} $Gui(WLA)
    pack $f.lNg -side top

    entry $f.eNg -width 5 -textvariable MarcacionElipsoide(Ng)
    bind $f.eNg <Return> "MarcacionElipsoideEstablece 5"
    pack $f.eNg -side top

    eval {label $f.lNs -text "Ns, FactConv y Carga de fichero"} $Gui(WLA)
    pack $f.lNs -side top
    
    frame $f.fNs -width 10
    entry $f.fNs.eNs -width 3 -textvariable MarcacionElipsoide(Ns)
    bind $f.fNs.eNs <Return> "MarcacionElipsoideEstablece 6"
    pack $f.fNs.eNs -side left -in $f.fNs
    entry $f.fNs.eFc -width 3 -textvariable MarcacionElipsoide(FactConv)
    bind $f.fNs.eFc <Return> "MarcacionElipsoideEstablece 26"
    pack $f.fNs.eFc -side left -in $f.fNs
    entry $f.fNs.eCa -width 3 -textvariable MarcacionElipsoide(Cargar)
    #bind $f.fNs.eFc <Return> "MarcacionElipsoideEstablece 26"
    pack $f.fNs.eCa -side left -in $f.fNs
    pack $f.fNs -side top

    eval {label $f.lVe -text "Valores del supervector de pesos"} $Gui(WLA)
    pack $f.lVe -side top

    frame $f.fVe -width 25
    entry $f.fVe.e0 -width 2 -textvariable MarcacionElipsoide(Ve0)
    bind $f.fVe.e0 <Return> "MarcacionElipsoideEstablece 7"
    pack $f.fVe.e0 -side left -in $f.fVe
    entry $f.fVe.e1 -width 2 -textvariable MarcacionElipsoide(Ve1)
    bind $f.fVe.e1 <Return> "MarcacionElipsoideEstablece 8"
    pack $f.fVe.e1 -side left -in $f.fVe
    entry $f.fVe.e2 -width 2 -textvariable MarcacionElipsoide(Ve2)
    bind $f.fVe.e2 <Return> "MarcacionElipsoideEstablece 9"
    pack $f.fVe.e2 -side left -in $f.fVe
    entry $f.fVe.e3 -width 2 -textvariable MarcacionElipsoide(Ve3)
    bind $f.fVe.e3 <Return> "MarcacionElipsoideEstablece 10"
    pack $f.fVe.e3 -side left -in $f.fVe
    entry $f.fVe.e4 -width 2 -textvariable MarcacionElipsoide(Ve4)
    bind $f.fVe.e4 <Return> "MarcacionElipsoideEstablece 11"
    pack $f.fVe.e4 -side left -in $f.fVe
    entry $f.fVe.e5 -width 2 -textvariable MarcacionElipsoide(Ve5)
    bind $f.fVe.e5 <Return> "MarcacionElipsoideEstablece 12"
    pack $f.fVe.e5 -side left -in $f.fVe
    entry $f.fVe.e6 -width 2 -textvariable MarcacionElipsoide(Ve6)
    bind $f.fVe.e6 <Return> "MarcacionElipsoideEstablece 13"
    pack $f.fVe.e6 -side left -in $f.fVe
    pack $f.fVe -side top

    eval {label $f.lPot -text "Valores de parámetros de funciones potenciales"} $Gui(WLA)
    pack $f.lPot -side top

    frame $f.fAlfa -width 25
    entry $f.fAlfa.ePot0 -width 2 -textvariable MarcacionElipsoide(AlfaPot0)
    bind $f.fAlfa.ePot0 <Return> "MarcacionElipsoideEstablece 14"
    pack $f.fAlfa.ePot0 -side left -in $f.fAlfa
    entry $f.fAlfa.ePot1 -width 2 -textvariable MarcacionElipsoide(AlfaPot1)
    bind $f.fAlfa.ePot1 <Return> "MarcacionElipsoideEstablece 15"
    pack $f.fAlfa.ePot1 -side left -in $f.fAlfa
    entry $f.fAlfa.ePot2 -width 2 -textvariable MarcacionElipsoide(AlfaPot2)
    bind $f.fAlfa.ePot2 <Return> "MarcacionElipsoideEstablece 16"
    pack $f.fAlfa.ePot2 -side left -in $f.fAlfa
    entry $f.fAlfa.ePot3 -width 2 -textvariable MarcacionElipsoide(AlfaPot3)
    bind $f.fAlfa.ePot3 <Return> "MarcacionElipsoideEstablece 17"
    pack $f.fAlfa.ePot3 -side left -in $f.fAlfa
    entry $f.fAlfa.ePot4 -width 2 -textvariable MarcacionElipsoide(AlfaPot4)
    bind $f.fAlfa.ePot4 <Return> "MarcacionElipsoideEstablece 18"
    pack $f.fAlfa.ePot4 -side left -in $f.fAlfa
    entry $f.fAlfa.ePot5 -width 3 -textvariable MarcacionElipsoide(AlfaPot5)
    bind $f.fAlfa.ePot5 <Return> "MarcacionElipsoideEstablece 19"
    pack $f.fAlfa.ePot5 -side left -in $f.fAlfa
    entry $f.fAlfa.ePot6 -width 2 -textvariable MarcacionElipsoide(AlfaPot6)
    bind $f.fAlfa.ePot6 <Return> "MarcacionElipsoideEstablece 20"
    pack $f.fAlfa.ePot6 -side left -in $f.fAlfa
    pack $f.fAlfa -side top

    frame $f.fBeta -width 23
    entry $f.fBeta.ePot0 -width 3 -textvariable MarcacionElipsoide(BetaPot0)
    bind $f.fBeta.ePot0 <Return> "MarcacionElipsoideEstablece 21"
    pack $f.fBeta.ePot0 -side left -in $f.fBeta
    entry $f.fBeta.ePot1 -width 3 -textvariable MarcacionElipsoide(BetaPot1)
    bind $f.fBeta.ePot1 <Return> "MarcacionElipsoideEstablece 22"
    pack $f.fBeta.ePot1 -side left -in $f.fBeta
    entry $f.fBeta.ePot2 -width 3 -textvariable MarcacionElipsoide(BetaPot2)
    bind $f.fBeta.ePot2 <Return> "MarcacionElipsoideEstablece 23"
    pack $f.fBeta.ePot2 -side left -in $f.fBeta
    entry $f.fBeta.ePot3 -width 3 -textvariable MarcacionElipsoide(BetaPot3)
    bind $f.fBeta.ePot3 <Return> "MarcacionElipsoideEstablece 24"
    pack $f.fBeta.ePot3 -side left -in $f.fBeta
    entry $f.fBeta.ePot4 -width 3 -textvariable MarcacionElipsoide(BetaPot4)
    bind $f.fBeta.ePot4 <Return> "MarcacionElipsoideEstablece 25"
    pack $f.fBeta.ePot4 -side left -in $f.fBeta
    pack $f.fBeta -side top
   
    
    eval {label $f.lcolor -text "MarcacionElipsoideEstablece color"} $Gui(WLA)
    pack $f.lcolor -side top

    entry $f.ecolor -width 4 -textvariable MarcacionElipsoide(color)
    pack $f.ecolor -side top

    #Añadimos un botón para calcular los contornos
    DevAddButton $f.bCalc "Calcular" "Tab MarcacionElipsoide row1 Calcular"
    TooltipAdd $f.bCalc "Go go go"
    pack $f.bCalc -side bottom -padx $Gui(pad) -pady $Gui(pad)
    
    #Añadimos un botón para trazar el elipsoide
    DevAddButton $f.bTraz "Trazo del elipsoide" "MarcacionElipsoideTrazar"
    pack $f.bTraz -side bottom

}
# end MarcacionElipsoideBuildParametrosFrame

#-------------------------------------------------------------------------------
# .PROC MarcacionElipsoideEstablece
#
# .END
#-------------------------------------------------------------------------------
proc MarcacionElipsoideEstablece {a} {

    global MarcacionElipsoide Module
    
	set id $MarcacionElipsoide(identif)

	if { $a==0 } then { MarcacionElipsoide(Generador,$id) SetTheta $MarcacionElipsoide(ResTheta) }
	if { $a==1 } then { MarcacionElipsoide(Generador,$id) SetPhi $MarcacionElipsoide(ResPhi) }
	if { $a==2 } then { MarcacionElipsoide(Plantilla,$id) SetJ $MarcacionElipsoide(J)
			       MarcacionElipsoide(Verosimil,$id) SetJ $MarcacionElipsoide(J) }
	if { $a==3 } then { MarcacionElipsoide(Plantilla,$id) SetK $MarcacionElipsoide(K)
			       MarcacionElipsoide(Verosimil,$id) SetK $MarcacionElipsoide(K)
			       MarcacionElipsoide(Optimiza,$id) SetK $MarcacionElipsoide(K) }
	if { $a==4 } then { MarcacionElipsoide(Plantilla,$id) Setdrmax $MarcacionElipsoide(drmax)
			       MarcacionElipsoide(Verosimil,$id) Setdrmax $MarcacionElipsoide(drmax) }
	if { $a==5 } then { MarcacionElipsoide(Verosimil,$id) SetNg $MarcacionElipsoide(Ng) }
	if { $a==27 } then { MarcacionElipsoide(Voluson,$id) SetFileName $MarcacionElipsoide(Nombre) 
			     MarcacionElipsoide(Plantilla,$id) SetFichModelo $MarcacionElipsoide(Nombre)
			     MarcacionElipsoide(Generador,$id) SetNomFichMar $MarcacionElipsoide(Nombre)
			     MarcacionElipsoide(Optimiza,$id) SetRuta $MarcacionElipsoide(Nombre) 
			     MarcacionElipsoide(Verosimil,$id) SetRuta $MarcacionElipsoide(Nombre) }
}
# end MarcacionElipsoideEstableceTheta


#-------------------------------------------------------------------------------
# .PROC MarcacionElipsoideBuildCalcularFrame
#
#   Crea el marco Calcular
# .END
#-------------------------------------------------------------------------------
proc MarcacionElipsoideBuildCalcularFrame {} {

    global Gui MarcacionElipsoide Module Volume

    #-------------------------------------------
    # Marco Calcular
    #-------------------------------------------
    set fCalcular $Module(MarcacionElipsoide,fCalcular)
    set f $fCalcular

    #Añadimos un botón
    DevAddButton $f.b1 "Cálculo elipsoide" "MarcacionElipsoideCalcular"
    pack $f.b1 -side top
    
    #Añadimos un botón
    DevAddButton $f.b2 "Cálculo modelo de verosimilitud" "MarcacionElipsoideModeloVerosimilitudCalc"
    pack $f.b2 -side top

    #Añadimos un botón
    DevAddButton $f.b3 "Cálculo optimización SA" "MarcacionElipsoideModeloOptimizacionSA"
    pack $f.b3 -side top
    
    #Añadimos un botón
    DevAddButton $f.b4 "Generar resultados validación" "ResultadosValidacion"
    pack $f.b4 -side top
    
    #Añadimos un botón
    DevAddButton $f.b5 "Dibujar resultado 2D" "DibujaResultado"
    pack $f.b5 -side top
    
    #Añadimos un botón
    DevAddButton $f.b6 "Dibujar un nuevo ajuste en 2D" "DibujaUnModelo"
    pack $f.b6 -side top
        
    #Añadimos un botón
    DevAddButton $f.b7 "Borrar todos los ajustes" "Borrado"
    pack $f.b7 -side top
    
    #Añadimos un botón
    DevAddButton $f.b8 "Regenerar todos los ajustes" "RegeneraTodos"
    pack $f.b8 -side top
    
    #Añadimos un botón
    DevAddButton $f.b9 "Borrar un ajuste" "BorradoModelo"
    pack $f.b9 -side top

    #Añadimos un botón
    DevAddButton $f.b10 "Cálculo completo" "RunMarcacionElipsoide"
    pack $f.b10 -side bottom
    
    #Añadimos un botón
    DevAddButton $f.b11 "Carga datos Voluson" "MarcacionElipsoideCarga"
    pack $f.b11 -side bottom

    #Añadimos un botón
    DevAddButton $f.b12 "Dibuja modelo en fichero" "MarcacionElipsoideDibModFich"
    pack $f.b12 -side bottom

}
# end MarcacionElipsoideBuildCalcularFrame


#-------------------------------------------------------------------------------
# .PROC MarcacionElipsoideEnter
# Llamado cuando el usuario entra en este módulo.  Pushes el manejador de eventos
# para este módulo. 
# .ARGS
# .END
#-------------------------------------------------------------------------------
proc MarcacionElipsoideEnter {} {
    global MarcacionElipsoide
    
    # Push el manejador de eventos
    #------------------------------------
    # Descripción:
    #   Para que los enlaces de eventos de este módulo no entren en conflicto con otros 
    #   módulos, se usan nuestros enlaces sólo cuando el usuario está en este módulo.
    #   La rutina pushEventManager salva los enlaces previos en 
    #   un stack y enlaza los nuevos nuestros.
    #   (Ver slicer/program/tcl-shared/Events.tcl para más detalles)
    pushEventManager $MarcacionElipsoide(eventManager)


    # limpiar el textBox e introducir las instrucciones allí
#    $MarcacionElipsoide(textBox) delete 1.0 end
#    $MarcacionElipsoide(textBox) insert end "Shift-Click anywhere!\n"

}


#-------------------------------------------------------------------------------
# .PROC MarcacionElipsoideExit
# Llamado cuando este el usuario sale del módulo. Pops el manejador de eventos
# para este módulo.  
# .ARGS
# .END
#-------------------------------------------------------------------------------
proc MarcacionElipsoideExit {} {

    # Pop el manejador de eventos
    #------------------------------------
    # Descripción:
    #   Usar esto con pushEventManager.  popEventManager borran nuestros 
    #   enlaces cuando el usuario sale del módulo y reubica 
    #   los anteriores.
    #
    popEventManager
}

#-------------------------------------------------------------------------------
# .PROC MarcacionElipsoideSetOriginal
# 
# .ARGS
# .END
#-------------------------------------------------------------------------------
proc MarcacionElipsoideSetOriginal {v} {
    global MarcacionElipsoide Volume
    
    if {$v == $MarcacionElipsoide(idWorking)} {
        tk_messageBox -message "The Original and Working volumes must differ."
        return
    }
    set MarcacionElipsoide(idOriginal) $v
    
    # Change button text
    $MarcacionElipsoide(mbOriginal) config -text [Volume($v,node) GetName]
}

#-------------------------------------------------------------------------------
# .PROC MarcacionElipsoideSetWorking
# 
# .ARGS
# .END
#-------------------------------------------------------------------------------
proc MarcacionElipsoideSetWorking {v} {
    global MarcacionElipsoide Volume Gui
    
    if {$v == [MarcacionElipsoideGetOriginalID]} {
        tk_messageBox -message "The Original and Working volumes must differ."
        return
    }

    set MarcacionElipsoide(idWorking) $v
    
    # Change button text, show name and file prefix
    if {$v == "NEW"} {
        $MarcacionElipsoide(mbWorking) config -text $v
        set MarcacionElipsoide(prefixWorking) ""
        set MarcacionElipsoide(nameWorking) Working
        eval {$MarcacionElipsoide(eNameWorking) configure -state normal}  $Gui(WEA)
    } else {
        $MarcacionElipsoide(mbWorking) config -text [Volume($v,node) GetName]
        set MarcacionElipsoide(prefixWorking) [MainFileGetRelativePrefix \
            [Volume($v,node) GetFilePrefix]]
        set MarcacionElipsoide(nameWorking) [Volume($v,node) GetName]
        # Disable name entry field if not NEW volume
        eval {$MarcacionElipsoide(eNameWorking) configure -state disabled} $Gui(WEDA)
    }

}

#-------------------------------------------------------------------------------
# .PROC MarcacionElipsoideGetOriginalID
# 
# .ARGS
# .END
#-------------------------------------------------------------------------------
proc MarcacionElipsoideGetOriginalID {} {
    global MarcacionElipsoide
    
    return $MarcacionElipsoide(idOriginal)
}

#-------------------------------------------------------------------------------
# .PROC MarcacionElipsoideGetWorkingID
#
# Returns the working volume's ID.
# If there is no working volume (MarcacionElipsoide(idWorking)==NEW), then it creates one.
# .END
#-------------------------------------------------------------------------------
proc MarcacionElipsoideGetWorkingID {} {
    global MarcacionElipsoide Volume Lut

    # If there is no Working volume, then create one
    if {$MarcacionElipsoide(idWorking) != "NEW"} {
        return $MarcacionElipsoide(idWorking)
    }
    
    # Create the node
    set n [MainMrmlAddNode Volume]
    set v [$n GetID]

    $n SetDescription "Working Volume=$v"
    $n SetLUTName     0
    $n InterpolateOff
    $n LabelMapOff
    
    # Make sure the name entered is okay, else use default
    if {[ValidateName $MarcacionElipsoide(nameWorking)] == 0} {
        tk_messageBox -message "The Descriptive Name can consist of letters, digits, dashes, or underscores only. Using default name Working"
        $n SetName Working
    } else {
        $n SetName $MarcacionElipsoide(nameWorking)   
    }
    
    # Create the volume
    MainVolumesCreate $v
    Volume($v,vol) UseLabelIndirectLUTOff

    MarcacionElipsoideSetWorking $v

    # This updates all the buttons to say that the
    # Volume List has changed.
    MainUpdateMRML

    return $v
}


#-------------------------------------------------------------------------------
# .PROC MarcacionElipsoideCount
#
# Esta rutina demuestra como hacer callbacks de botones y usar arrays globales
# para programación orientada a objetos.
# .END
#-------------------------------------------------------------------------------
proc MarcacionElipsoideCount {} {
    global MarcacionElipsoide
    
    incr MarcacionElipsoide(count)
    $MarcacionElipsoide(lParameters) config -text "Has pulsado el botón $MarcacionElipsoide(count) veces"
}


#-------------------------------------------------------------------------------
# .PROC MarcacionElipsoideShowFile
#
# Esta rutina demuestra como hacer callbacks de botones y usar arrays globales
# para programación orientada a objetos.
# .END
#-------------------------------------------------------------------------------
proc MarcacionElipsoideShowFile {} {
    global MarcacionElipsoide
    
    $MarcacionElipsoide(lfile) config -text "Has entrado en: $MarcacionElipsoide(FileName)"
}


#-------------------------------------------------------------------------------
# .PROC MarcacionElipsoideBindingCallback
# Demostración de rutinas callback para enlaces
# 
# .ARGS
# .END
#-------------------------------------------------------------------------------
proc MarcacionElipsoideBindingCallback { event W X Y x y t } {
    global MarcacionElipsoide

    set insertText "$evento en: $X $Y\n"
    
    switch $event {
    "Shift-2" {
        set insertText "Don't poke the Slicer!\n"
    }
    "Shift-3" {
        set insertText "Ouch!\n"
    }

    }
#    $MarcacionElipsoide(textBox) insert end $insertText

}


#-------------------------------------------------------------------------------
# .PROC MarcacionElipsoidePrepareResult
#   Crea el nuevo volumen si es necesario. En otro caso, pregunta para sobreescribir.
#   Devuelve 1 si hay errores y 0 en otro caso.
# .END
#-------------------------------------------------------------------------------
proc MarcacionElipsoideCheckErrors {} {
    global MarcacionElipsoide Volume

    if {  ($MarcacionElipsoide(InputVol) == $Volume(idNone)) || \
          ($MarcacionElipsoide(ResultVol)   == $Volume(idNone))}  {
        DevErrorWindow "No se puede usar el Volumen \"None\""
        return 1
    }
    return 0
}


#-------------------------------------------------------------------------------
# .PROC MarcacionElipsoidePrepareResultVolume
#   Chequea errores en la inicialización.
#   Devuelve 1 si hay errores y 0 en otro caso.
# .END
#-------------------------------------------------------------------------------
proc MarcacionElipsoidePrepareResultVolume {}  {
    global MarcacionElipsoide

    set v1 $MarcacionElipsoide(InputVol)
    set v2 $MarcacionElipsoide(ResultVol)

    # Necesitamos crear un nuevo volumen?
    # Si es así, hagámoslo.
    
    if {$v2 == -5 } {
        set v2 [DevCreateNewCopiedVolume $v1 ""  "MarcacionElipsoideResult" ]
        set node [Volume($v2,vol) GetMrmlNode]
        Mrml(dataTree) RemoveItem $node 
        set nodeBefore [Volume($v1,vol) GetMrmlNode]
    Mrml(dataTree) InsertAfterItem $nodeBefore $node
        MainUpdateMRML
    } else {

        # Estamos sobreescribiendo un volumen?
        # Si es así, pregúntese; si no, sálgase.
         
        set v2name  [Volume($v2,node) GetName]
    set continue [DevOKCancel "Overwrite $v2name?"]
          
        if {$continue == "cancel"} { return 1 }
        # ¡Han dicho que OK, así que sobreescribamos!
              
        Volume($v2,node) Copy Volume($v1,node)
    }

    set MarcacionElipsoide(ResultVol) $v2
    

    return 0
}


#-------------------------------------------------------------------------------
# .PROC RunMarcacionElipsoide
#   Ejecutar el algoritmo completo
#
# .END
#-------------------------------------------------------------------------------
proc RunMarcacionElipsoide {} {

  global MarcacionElipsoide Module Volume Gui

  MarcacionElipsoideCalcular 
  MarcacionElipsoideModeloVerosimilitudCalc
  MarcacionElipsoideModeloOptimizacionSA

}

#-------------------------------------------------------------------------------
# .PROC MarcacionElipsoideDibModFich
#   Lee el modelo de un fichero y lo dibuja
#
# .END
#-------------------------------------------------------------------------------

proc MarcacionElipsoideDibModFich {} {

global MarcacionElipsoide Module Volume Gui View

	set id $MarcacionElipsoide(identif)
	set tipo 0
	set MarcacionElipsoide(color) 0
	
	set MarcacionElipsoide(Poli,$MarcacionElipsoide(color)) [MarcacionElipsoide(Plantilla,$id) LeeModelo]

	set max [MarcacionElipsoide(Plantilla,$id) GetMax]
	set min [MarcacionElipsoide(Plantilla,$id) GetMin]
	MarcacionElipsoide(Generador,$id) SetMax [lindex $max 0] [lindex $max 1] [lindex $max 2] [lindex $max 3]
	MarcacionElipsoide(Generador,$id) SetMin [lindex $min 0] [lindex $min 1] [lindex $min 2] [lindex $min 3]

	set MarcacionElipsoide(actor,$MarcacionElipsoide(color)) [MarcacionElipsoide(Generador,$id) DibujaModelo $MarcacionElipsoide(Poli,$MarcacionElipsoide(color)) $tipo $MarcacionElipsoide(color)]	
	MainAddActor $MarcacionElipsoide(actor,$MarcacionElipsoide(color))

	if { $tipo==1 } { 
	set MarcacionElipsoide(actor,[expr 1+$MarcacionElipsoide(color)]) [MarcacionElipsoide(Generador,$id) DevuelveBarra]
	MainAddActor $MarcacionElipsoide(actor,[expr 1+$MarcacionElipsoide(color)]) }
	
  	RenderAll
}

#-------------------------------------------------------------------------------
# .PROC RunMarcacionElipsoide
#   Ejecutar el algoritmo completo
#
# .END
#-------------------------------------------------------------------------------

proc MarcacionElipsoideCarga {} {

global MarcacionElipsoide Module Volume Gui View

    set id $MarcacionElipsoide(identif)

    set w [MarcacionElipsoideGetWorkingID]
	  set prov [MarcacionElipsoide(Voluson,$id) GetOutput]
	  $prov Update
	  puts [$prov GetSpacing]
	  puts [$prov GetDimensions]
	  eval Volume($w,node) SetSpacing [$prov GetSpacing]
	  eval Volume($w,node) SetDimensions [lindex [$prov GetDimensions] 0] [lindex [$prov GetDimensions] 1]
	  Volume($w,node) SetImageRange 1 [lindex [$prov GetDimensions] 2]
	  Volume($w,vol) SetImageData $prov
	
	  set fov 0
  	for {set i 0} {$i < 2} {incr i} {
	    set dim     [lindex [Volume($w,node) GetDimensions] $i]
       	    set spacing [lindex [Volume($w,node) GetSpacing] $i]
           set newfov     [expr $dim * $spacing]
       if { $newfov > $fov } {
            set fov $newfov
            puts $fov
       }
    }
    set dim [lindex [Volume($w,node) GetImageRange] 1]
    set spacing [lindex [Volume($w,node) GetSpacing] 2]
    set newfov [expr $dim * $spacing]
    if { $newfov > $fov } {
          set fov $newfov 
    }
    puts $fov
  	set View(fov) $fov
		if { $MarcacionElipsoide(CambiaFOV)==1 } { 
  	MainViewSetFov
    set MarcacionElipsoide(CambiaFOV) 0 }
  	Volume($w,node) ComputeRasToIjkFromScanOrder $Volume(scanOrder)
    MainVolumesUpdate $w
    MainUpdateMRML
    MainSlicesSetVolumeAll Back $w
    #SegmentaMiocardio(Generador,$id) IntroduceMatriz [Volume($w,node) GetRasToIjk]
    
}

#-------------------------------------------------------------------------------
# .PROC MarcacionElipsoideTrazar
#   Traza el elipsoide para visualización de que son correctos los puntos tomados
#
# .END
#-------------------------------------------------------------------------------
proc MarcacionElipsoideTrazar {} {

  global MarcacionElipsoide Module Volume Gui

	set id $MarcacionElipsoide(identif)
	set poli [MarcacionElipsoide(Generador,$id) GeneraElipsoide Fiducials($id,points)]
	set color 0
	set tipo 0
	set MarcacionElipsoide(actora) [MarcacionElipsoide(Generador,$id) DibujaModelo $poli $tipo $color]	
	MainAddActor $MarcacionElipsoide(actora)
	if { $tipo==1 } { 
	set MarcacionElipsoide(actorb) [MarcacionElipsoide(Generador,$id) DevuelveBarra]
	MainAddActor $MarcacionElipsoide(actorb) }

	
	Render3D	
}

#-------------------------------------------------------------------------------
# .PROC DibujaResultado
#   Dibuja el resultado y el ajuste original en 2D
#
# .END
#-------------------------------------------------------------------------------
proc DibujaResultado {} {

  global MarcacionElipsoide Module Volume Gui Ed
  	
	set id $MarcacionElipsoide(identif)
			
	set plan [Slicer GetReformatMatrix 1]
	set e [$plan GetElement 0 3]
	
	set piwi [lindex [Volume($MarcacionElipsoide(idWorking),node) GetSpacing] 0]
	set pihe [lindex [Volume($MarcacionElipsoide(idWorking),node) GetSpacing] 1]
	#set slth [lindex [Volume($MarcacionElipsoide(idWorking),node) GetSpacing] 2]
	#tk_messageBox -message "$piwi"
	#tk_messageBox -message "$pihe"

	set M [MarcacionElipsoide(Verosimil,$id) GetMbis]
	set N [MarcacionElipsoide(Verosimil,$id) GetNbis]
	if { $N > $M } then { set M $N }
	#set matriz [Volume($MarcacionElipsoide(idWorking),node) GetRasToIjk]
	set centro [MarcacionElipsoide(Plantilla,$id) ObtieneCentro]
	set puntoconcretoa [$centro GetComponent 1 0]
	set puntoconcretob [$centro GetComponent 0 0]
	set incrz [expr $puntoconcretoa-$puntoconcretob]
	set z [expr int(($e-$puntoconcretob)/$incrz)]
	#tk_messageBox -message "$z"

	Slicer DrawDeleteAll
	Slicer SetActiveSlice 1
	Slicer DrawSetShapeToLines
	Slicer DrawShowPoints 1
	Slicer DrawSetColor 1 0 0
	Slicer DrawSetRadius 0
	
	for { set a 0 } { $a < $MarcacionElipsoide(J) } { incr a } {
		
		set x [lindex [$MarcacionElipsoide(Poli,0) GetPoint [expr $z*$MarcacionElipsoide(J)+$a]] 1]
		set y [lindex [$MarcacionElipsoide(Poli,0) GetPoint [expr $z*$MarcacionElipsoide(J)+$a]] 2]
		#tk_messageBox -message "$x"
		Slicer DrawInsertPoint [expr int((-$x/$piwi+$M/2)*256/$M)] [expr int(($y/$pihe+$M/2)*256/$M)] }
	
	Slicer DrawComputeIjkPoints
    	set points [Slicer GetDrawIjkPoints]
    	set v [MarcacionElipsoideGetWorkingID]
    	MarcacionElipsoideSetupBeforeApplyEffect $v Single Active
    	Ed(editor) Draw 5 $points 0 Lines
	Ed(editor)   SetInput ""
    	Ed(editor)   UseInputOff
	Slicer DrawDeleteAll
	#catch {__EditorPending_Points Delete}
	MarcacionElipsoideUpdateAfterApplyEffect $v Active
	
	Slicer DrawSetShapeToLines
	Slicer DrawShowPoints 1
	Slicer DrawSetColor 0 1 0
	Slicer DrawSetRadius 0
	
	for { set b 0 } { $b < $MarcacionElipsoide(J) } { incr b } {
		set x [lindex [$MarcacionElipsoide(Poli,1) GetPoint [expr $z*$MarcacionElipsoide(J)+$b]] 1]
		set y [lindex [$MarcacionElipsoide(Poli,1) GetPoint [expr $z*$MarcacionElipsoide(J)+$b]] 2]
		Slicer DrawInsertPoint [expr int((-$x/$piwi+$M/2)*256/$M)] [expr int(($y/$pihe+$M/2)*256/$M)] }
		
	Slicer DrawComputeIjkPoints
    	set pointsbis [Slicer GetDrawIjkPoints]
    	MarcacionElipsoideSetupBeforeApplyEffect $v Single Active
    	Ed(editor) Draw 6 $pointsbis 0 Lines
    	Ed(editor)   SetInput ""
    	Ed(editor)   UseInputOff
	Slicer DrawDeleteAll
	#catch {__EditorPending_Points Delete}
	MarcacionElipsoideUpdateAfterApplyEffect $v Active
		
}


#-------------------------------------------------------------------------------
# .PROC DibujaUnModelo
#   Dibuja un modelo (dado por el color) en 2D
#
# .END
#-------------------------------------------------------------------------------
proc DibujaUnModelo {} {

  global MarcacionElipsoide Module Volume Gui Ed
  	
	set id $MarcacionElipsoide(identif)
			
	set plan [Slicer GetReformatMatrix 1]
	set e [$plan GetElement 0 3]
	
	set piwi [lindex [Volume($MarcacionElipsoide(idWorking),node) GetSpacing] 0]
	set pihe [lindex [Volume($MarcacionElipsoide(idWorking),node) GetSpacing] 1]
	#set slth [lindex [Volume($MarcacionElipsoide(idWorking),node) GetSpacing] 2]
	#tk_messageBox -message "$piwi"
	#tk_messageBox -message "$pihe"

	set M [MarcacionElipsoide(Verosimil,$id) GetMbis]
	set N [MarcacionElipsoide(Verosimil,$id) GetNbis]
	if { $N > $M } then { set M $N }
	#set matriz [Volume($MarcacionElipsoide(idWorking),node) GetRasToIjk]
	set centro [MarcacionElipsoide(Plantilla,$id) ObtieneCentro]
	set puntoconcretoa [$centro GetComponent 1 0]
	set puntoconcretob [$centro GetComponent 0 0]
	set incrz [expr $puntoconcretoa-$puntoconcretob]
	set z [expr int(($e-$puntoconcretob)/$incrz)]
	#tk_messageBox -message "$z"

	Slicer DrawDeleteAll
	Slicer SetActiveSlice 1
	Slicer DrawSetShapeToLines
	Slicer DrawShowPoints 1
	Slicer DrawSetColor 0 0 1
	Slicer DrawSetRadius 0
	
	for { set a 0 } { $a < $MarcacionElipsoide(J) } { incr a } {
		
		set x [lindex [$MarcacionElipsoide(Poli,$MarcacionElipsoide(color)) GetPoint [expr $z*$MarcacionElipsoide(J)+$a]] 1]
		set y [lindex [$MarcacionElipsoide(Poli,$MarcacionElipsoide(color)) GetPoint [expr $z*$MarcacionElipsoide(J)+$a]] 2]
		#tk_messageBox -message "$x"
		Slicer DrawInsertPoint [expr int((-$x/$piwi+$M/2)*256/$M)] [expr int(($y/$pihe+$M/2)*256/$M)] }
	
	#Slicer DrawComputeIjkPoints
    	#set points [Slicer GetDrawIjkPoints]
    	#set v [MarcacionElipsoideGetWorkingID]
    	#MarcacionElipsoideSetupBeforeApplyEffect $v Single Active
    	#Ed(editor) Draw 10 $points 0 Lines
	#Ed(editor)   SetInput ""
    	#Ed(editor)   UseInputOff
	#Slicer DrawDeleteAll
	#catch {__EditorPending_Points Delete}
	#MarcacionElipsoideUpdateAfterApplyEffect $v Active
}


#-------------------------------------------------------------------------------
# .PROC MarcacionElipsoideCalcular
#   Calcula el ajuste inicial para cada slice y rayo y lo posiciona
#
# .END
#-------------------------------------------------------------------------------
proc MarcacionElipsoideCalcular {} {

  global MarcacionElipsoide Module Volume Gui

	set id $MarcacionElipsoide(identif)

	if {$MarcacionElipsoide(Cargar) == 1} {
        set poli [MarcacionElipsoide(Generador,$id) GeneraElipsoide [MarcacionElipsoide(Generador,$id) CargarMarcacionInicial] ]
      } else {
	  set poli [MarcacionElipsoide(Generador,$id) GeneraElipsoide Fiducials($id,points)]
	  MarcacionElipsoide(Generador,$id) GuardarMarcacionInicial Fiducials($id,points)
      }

	set plan [Slicer GetReformatMatrix 1]

	MarcacionElipsoide(Plantilla,$id) IntroducePlano $plan
	set slth [lindex [Volume($MarcacionElipsoide(idWorking),node) GetSpacing] 2]
	MarcacionElipsoide(Plantilla,$id) SetSlth $slth
	MarcacionElipsoide(Plantilla,$id) Setdrmax [expr 18.75*$slth]
	
	MarcacionElipsoide(Plantilla,$id) ObtieneParametros $poli
}


#-------------------------------------------------------------------------------
# .PROC MarcacionElipsoideModeloVerosimilitudCalc
#   Calcula los términos energéticos del Campo de Gibbs correspondientes al modelo de 
#   verosimilitud
#
# .END
#-------------------------------------------------------------------------------
proc MarcacionElipsoideModeloVerosimilitudCalc {} {

  global MarcacionElipsoide Module Volume Gui

	set id $MarcacionElipsoide(identif)

      set id $MarcacionElipsoide(identif)
    
      MarcacionElipsoide(Verosimil,$id) SetP [MarcacionElipsoide(Plantilla,$id) GetP]
      set piwi [lindex [Volume($MarcacionElipsoide(idWorking),node) GetSpacing] 0]
      set pihe [lindex [Volume($MarcacionElipsoide(idWorking),node) GetSpacing] 1]
      set slth [lindex [Volume($MarcacionElipsoide(idWorking),node) GetSpacing] 2]
      MarcacionElipsoide(Verosimil,$id) Setdrmax [expr 18.75*$slth]
      set plan [Slicer GetReformatMatrix 1]
      MarcacionElipsoide(Verosimil,$id) IntroducePlano $plan
      MarcacionElipsoide(Verosimil,$id) IntroduceMatriza [Volume($MarcacionElipsoide(idWorking),node) GetRasToIjk]
      set lim [MarcacionElipsoide(Plantilla,$id) Getlimite]
      MarcacionElipsoide(Verosimil,$id) Setlimite [lindex $lim 0] [lindex $lim 1] [lindex $lim 2] [lindex $lim 3] [lindex $lim 4] [lindex $lim 5]
	MarcacionElipsoide(Verosimil,$id) IntroduceImagen [Volume($MarcacionElipsoide(idWorking),vol) GetOutput] 0
	MarcacionElipsoide(Verosimil,$id) IntroduceParam [MarcacionElipsoide(Plantilla,$id) ObtieneParam]
	MarcacionElipsoide(Verosimil,$id) IntroduceCentro [MarcacionElipsoide(Plantilla,$id) ObtieneCentro]
	MarcacionElipsoide(Verosimil,$id) Ejecuta
	#MarcacionElipsoide(Verosimil,$id) EscribeImagen
	
	MarcacionElipsoide(Plantilla,$id) EstableceMascaraElipsoide [MarcacionElipsoide(Verosimil,$id) DevuelveMascaraElipsoide]
	
	set tipo 0
	set color 0
	
	set rho [MarcacionElipsoide(Plantilla,$id) GeneraRhoNulo]
	set MarcacionElipsoide(Poli,$color) [MarcacionElipsoide(Plantilla,$id) ConstruyeModelo $rho]

	set MarcacionElipsoide(actor,$color) [MarcacionElipsoide(Generador,$id) DibujaModelo $MarcacionElipsoide(Poli,$color) $tipo $color]
	MainAddActor $MarcacionElipsoide(actor,$color)
	if { $tipo==1 } { 
	set MarcacionElipsoide(actor,[expr $color+1]) [MarcacionElipsoide(Generador,$id) DevuelveBarra]
	MainAddActor $MarcacionElipsoide(actor,[expr $color+1]) }

	RenderAll
}


#-------------------------------------------------------------------------------
# .PROC MarcacionElipsoideModeloOptimizacionSA
#   Calcula la optimización del contorno según el algoritmo SA
#
# .END
#-------------------------------------------------------------------------------
proc MarcacionElipsoideModeloOptimizacionSA {} {

  global MarcacionElipsoide Module Volume Gui
  
        set id $MarcacionElipsoide(identif)

	MarcacionElipsoide(Optimiza,$id) SetSuper 0
	MarcacionElipsoide(Optimiza,$id) SetPeriodMalla 1 0 0
	MarcacionElipsoide(Optimiza,$id) SetNumEntidades 1
	MarcacionElipsoide(Optimiza,$id) SetDimensionalidadEstados 1
	MarcacionElipsoide(Optimiza,$id) SetDimensionalidadMalla 2
	MarcacionElipsoide(Optimiza,$id) SetIndependencia 0
	MarcacionElipsoide(Optimiza,$id) SetOrdenMalla 4
	MarcacionElipsoide(Optimiza,$id) SetOrdenEstados 5

	MarcacionElipsoide(Optimiza,$id) SetDimMalla [MarcacionElipsoide(Verosimil,$id) GetJ] [MarcacionElipsoide(Verosimil,$id) GetP] 1
	MarcacionElipsoide(Optimiza,$id) SetDimZ 1 1 1
	MarcacionElipsoide(Optimiza,$id) SetK [MarcacionElipsoide(Verosimil,$id) GetK]
	MarcacionElipsoide(Optimiza,$id) SetV 2 
	MarcacionElipsoide(Optimiza,$id) ConstruyeVecindario
	MarcacionElipsoide(Optimiza,$id) EstableceLE [MarcacionElipsoide(Verosimil,$id) DevuelveLR] 0 0

	for { set b 0 } { $b < 1 } { incr b } {
	MarcacionElipsoide(Optimiza,$id) Llamada $b
	set MarcacionElipsoide(Poli,$MarcacionElipsoide(color)) [MarcacionElipsoide(Plantilla,$id) ConstruyeModelo [MarcacionElipsoide(Optimiza,$id) DevuelveRho]]
	ResultadosValidacion }
	set max [MarcacionElipsoide(Plantilla,$id) GetMax]
	set min [MarcacionElipsoide(Plantilla,$id) GetMin]
	MarcacionElipsoide(Generador,$id) SetMax [lindex $max 0] [lindex $max 1] [lindex $max 2] [lindex $max 3]
	MarcacionElipsoide(Generador,$id) SetMin [lindex $min 0] [lindex $min 1] [lindex $min 2] [lindex $min 3]
	set tipo 0
	set MarcacionElipsoide(actor,$MarcacionElipsoide(color)) [MarcacionElipsoide(Generador,$id) DibujaModelo $MarcacionElipsoide(Poli,$MarcacionElipsoide(color)) $tipo $MarcacionElipsoide(color)]	
	MainAddActor $MarcacionElipsoide(actor,$MarcacionElipsoide(color))

	if { $tipo==1 } { 
	set MarcacionElipsoide(actor,[expr 1+$MarcacionElipsoide(color)]) [MarcacionElipsoide(Generador,$id) DevuelveBarra]
	MainAddActor $MarcacionElipsoide(actor,[expr 1+$MarcacionElipsoide(color)]) }
	
  	RenderAll
}

#-------------------------------------------------------------------------------
# .PROC Borrado
#   Borra los actores
#
# .END
#-------------------------------------------------------------------------------
proc Borrado {} {

  global MarcacionElipsoide Module Volume Gui
  
 	for { set i 0 } { $i< [expr $MarcacionElipsoide(color)+1] } { incr i } { MainRemoveActor $MarcacionElipsoide(actor,$i) }
	RenderAll
}

#-------------------------------------------------------------------------------
# .PROC RegeneraTodos
#
# .END
#-------------------------------------------------------------------------------
proc RegeneraTodos {} {

  global MarcacionElipsoide Module Volume Gui
  
 	for { set i 0 } { $i< [expr $MarcacionElipsoide(color)+1] } { incr i } { MainAddActor $MarcacionElipsoide(actor,$i) }
	RenderAll
}

#-------------------------------------------------------------------------------
# .PROC ResultadosValidacion
#
# .END
#-------------------------------------------------------------------------------
proc ResultadosValidacion {} {

  global MarcacionElipsoide Module Volume Gui
   
   set id $MarcacionElipsoide(identif)
   MarcacionElipsoide(Valida,$id) IntroducePlantilla $MarcacionElipsoide(Poli,$MarcacionElipsoide(color))
   MarcacionElipsoide(Valida,$id) IntroduceMatriza [Volume($MarcacionElipsoide(idWorking),node) GetRasToIjk]
   MarcacionElipsoide(Valida,$id) GeneraValores
   set ejes [MarcacionElipsoide(Valida,$id) Ejecuta]
   #set MarcacionElipsoide(actor,2) [MarcacionElipsoide(Generador,$id) DibujaModelo $ejes 2]
	
   #MainAddActor $MarcacionElipsoide(actor,2)
	
   #RenderAll

}

#-------------------------------------------------------------------------------
# .PROC BorradoModelo
#   Borra un actor correspondiente a un modelo (el dado por el color)
#
# .END
#-------------------------------------------------------------------------------
proc BorradoModelo {} {

  global MarcacionElipsoide Ed Module Volume Gui Editor 

 	MainRemoveActor $MarcacionElipsoide(actor,$MarcacionElipsoide(color))
 	unset MarcacionElipsoide(actor,$MarcacionElipsoide(color))
 	unset MarcacionElipsoide(Poli,$MarcacionElipsoide(color))
	RenderAll
}


#-------------------------------------------------------------------------------
# .PROC MarcacionElipsoideUpdateAfterApplyEffect
# 
# .ARGS
# .END
#-------------------------------------------------------------------------------
proc MarcacionElipsoideUpdateAfterApplyEffect {v {render All}} {
    global Ed Volume Lut Editor Slice
    
    set o [MarcacionElipsoideGetOriginalID]
    set w [MarcacionElipsoideGetWorkingID]
    
    # Get output from editor
    Volume($w,vol) SetImageData [Ed(editor) GetOutput]
    EditorActivateUndo [Ed(editor) GetUndoable]
    
    # w copies o's MrmlNode if the Input was the Original
    if {$v == $o} {
        EditorCopyNode $w $o
    }
    
    # Keep a copy for undo
    Editor(undoNode) Copy Volume($w,node)
    
    # Update pipeline and GUI
    MainVolumesUpdate $w
    
    # Render
    Render$render
    
    # Mark the volume as changed
    set Volume($w,dirty) 1
    
    $Editor(lRunTime)   config -text \
        "[format "%.2f" [Ed(editor) GetRunTime]] sec,"
    $Editor(lTotalTime) config -text \
        "[format "%.2f" [Ed(editor) GetTotalTime]] sec"
}

#-------------------------------------------------------------------------------
# .PROC MarcacionElipsoideSetupBeforeApplyEffect
# 
# .ARGS
# .END
#-------------------------------------------------------------------------------
proc MarcacionElipsoideSetupBeforeApplyEffect {v scope multi} {
    global Volume Ed Editor Gui
    
    set o [MarcacionElipsoideGetOriginalID]
    set w [MarcacionElipsoideGetWorkingID]
    
    if {[EditorSameExtents $w $o] != 1} {
        EditorCopyNode $w $o
        MainVolumesCopyData $w $o On
        # >> AT 12/20/01
        # force the working volume to be of type short
        set workvol [Volume($w,vol) GetOutput]
        set worktype [$workvol GetScalarType]
        # 4 is the ID of short in VTK
        if {$worktype != "4"} {
            $workvol SetScalarType 4
            $workvol AllocateScalars
        }
        # << AT 12/20/01
	#EditorClear Working        
    }
    
    # Set the editor's input & output
    Ed(editor) SetInput [Volume($o,vol) GetOutput]
    if {$v == $w} {
        Ed(editor) SetOutput [Volume($w,vol) GetOutput]
        Ed(editor) UseInputOff
    } else {
        Ed(editor) UseInputOn
    }
    
    Ed(editor) SetDimensionTo$scope
    
    # Set the slice orientation and number
    # (not used for 3D)
    
    set s      [Slicer GetActiveSlice]
    set orient [Slicer GetOrientString $s]
    set slice  [Slicer GetOffset $s]
    
    if {[lsearch "AxiSlice CorSlice SagSlice" $orient] == -1} {
        tk_messageBox -icon warning -title $Gui(title) -message \
            "The orientation of the active slice\n\
            must be one of: AxiSlice, CorSlice, SagSlice"
        return
    }
    switch $orient {
        "AxiSlice" {
            set order IS
        }
        "SagSlice" {
            set order LR
        }
        "CorSlice" {
            set order PA
        }
    }
    
    # Does the user want the orien of the active slice or native slices?
    if {$scope == "Multi" && $multi == "Native"} {
        set order [Volume($o,node) GetScanOrder]
    }
    switch $order {
        "SI" {
            set order IS
        }
        "RL" {
            set order LR
        }
        "AP" {
            set order PA
        }
    }
    
    Ed(editor) SetOutputSliceOrder $order
    Ed(editor) SetInputSliceOrder [Volume($v,node) GetScanOrder]
    Ed(editor) SetSlice $slice
}

#-------------------------------------------------------------------------------
# .PROC EditorClear
#
# Clear either the Working or Composite data to all zeros.
# .END
#-------------------------------------------------------------------------------
proc MarcacionElipsoideClear {data} {
    global Volume Editor Slice
    
    # If the volume doesn't exist yet, then don't write it, duh!
    if {$Editor(id$data) == "NEW"} {
        tk_messageBox -message "Nothing to clear."
        return
    }
    
    switch $data {
        Composite {set v $Editor(idComposite)}
        Working   {set v $Editor(idWorking)}
    }
    
    vtkImageCopy copy
    copy ClearOn
    copy SetInput [Volume($v,vol) GetOutput]
    copy Update
    copy SetInput ""
    Volume($v,vol) SetImageData [copy GetOutput]
    copy SetOutput ""
    copy Delete
    
    # Mark the volume as changed
    set Volume($v,dirty) 1
    
    MainVolumesUpdate $v
    RenderAll
}

#-------------------------------------------------------------------------------
# .PROC MarcacionElipsoideVolumenFiltrado
#
# Añade un volumen filtrado correspondiente al que se quiere segmentar (para imagen de borde).
# .END
#-------------------------------------------------------------------------------
proc MarcacionElipsoideVolumenFiltrado { } {
    
    global MarcacionElipsoide Volume Editor Slice
    
    set id $MarcacionElipsoide(identif)
    
    MarcacionElipsoide(Verosimil,$id) SetP [MarcacionElipsoide(Plantilla,$id) GetP]
    set piwi [lindex [Volume($MarcacionElipsoide(idWorking),node) GetSpacing] 0]
    set pihe [lindex [Volume($MarcacionElipsoide(idWorking),node) GetSpacing] 1]
    set slth [lindex [Volume($MarcacionElipsoide(idWorking),node) GetSpacing] 2]
    set plan [Slicer GetReformatMatrix 1]
    MarcacionElipsoide(Verosimil,$id) IntroducePlano $plan
    MarcacionElipsoide(Verosimil,$id) IntroduceMatriza [Volume($MarcacionElipsoide(idWorking),node) GetRasToIjk]
    MarcacionElipsoide(Verosimil,$id) Setlimite [MarcacionElipsoide(Plantilla,$id) Getlimite] 
    MarcacionElipsoide(Verosimil,$id) IntroduceImagen [Volume($MarcacionElipsoide(idWorking),vol) GetOutput] 1   
}

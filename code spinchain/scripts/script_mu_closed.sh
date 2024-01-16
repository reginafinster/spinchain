#!/bin/bash
CURRENTPATH="$PWD"
#CURRENTPATH='/home/rfinster/Masterarbeit/Code/code_aktueller_Stand_neu/closed_chain_transient/'


for realization in {1,2,3,4,5,6}
do
for disorder in {10,100}
do

RUNPATH="$CURRENTPATH/h_$disorder/realization_$realization"
mkdir -p $RUNPATH
param=`echo "scale=3; $disorder/100"|bc -l`
sed "s/#disorder/$param/g" < parameters > closed_chain_random_h${param}.cfg
mv closed_chain_random_h${param}.cfg $RUNPATH
cp heisenberg_liouville_space $RUNPATH
cd $RUNPATH
echo "realization $realization disorder $disorder"
# ./heisenberg_liouville_space open_chain_random_mu${param}.cfg
qsub -mem 2 -args closed_chain_random_h${param}.cfg heisenberg_liouville_space
sleep 1
cd $CURRENTPATH

done
done
done








# echo gibt auf shell aus
# das $-Zeichen gibt bei echo Variablen aus
# Rechnen in der Shell: echo $(( (23-2)*2/3 )), also in Klammern hinter echo kann eine einfache Rechenoperation ausgeführt werden, diese wird dann intern wohl als Variable gespeichert und benötigt daher das Dollarzeichen zur Ausgabe.
# Rechenbefehl einleiten mit doppelten runden Klammern
# aber die shell kann nur mit integern rechnen
# der pipe-operator | leitet die Ausgabe des einen Befehls an den nächsten Befehl weiter
# dieser <<< operator leitet einen string auf die stdin eines Befehls weiter.

# bc ist ein interaktives Rechenprogramm, das auch Fließkommazahlen kann. Die shell kann wirklich nur integer verarbeiten, dh. ich muss alle Schleifen dann mit integern aufsetzen und danach mit bc umrechnen.

# sed ist ein stream editor, mit dem man Text streamen und dabei bearbeiten kann (zum Beispiel kann man ersetzungen vornehmen). man ruft ihn auf mit
# sed sed-skript textdatei. Das Ergebnis wird per default auf die stdin ausgegeben. 
# s/#mu/$param/g bedeutet: jedes Auftreten von #mu wird durch den parameter param ersetzt
# sed s/Anton/Berta/g Textdatei ersetzt Anton durch Bertha in Textdatei.
# Ich denke hier wird der String verwendet weil zuerst innerhalb des Skripts die ersetzung für den parameter durchgeführt werden muss. Dann braucht man den Umleitungsoperator, der die Eingabe auf die Datei parameters umleitet 

# > dieser Befehl leitet z.b. in eine Datei um
# < dieser Befehl leitet die Standardeingabe um, z.b. tr -d '0-9' < datei.txt zeigt den Inhalt von datei.txt ohne ziffern an

# Kelvin-Voigt-Cardiomyocyte

Zweidimensionales Kelvin-Voigt-Modell mit Trägheit und aktiven
longitudinalen Sarkomerkräften, implementiert in R.

Das Projekt formuliert die Modellgleichungen zunächst symbolisch,
diskretisiert sie auf einem räumlichen Gitter und simuliert das
resultierende ODE-System mit dMod2.

## Modell

Die Zustandsgrößen sind die vektorwertige Verschiebung u und
die Geschwindigkeit v:

    du/dt = v

    rho * dv/dt = div(E : grad(u) + D : grad(v)) - kappa*u + P

Dabei stehen:

- E: Elastizitätstensor
- D: Viskositätstensor
- rho: Massendichte
- kappa: lokale Rückstellsteifigkeit
- P: aktive Kraftdichte

Die aktive Kraft wird aus paarweise entgegengesetzten Zugkräften
zwischen longitudinal benachbarten Knoten aufgebaut.

In der Beispielsimulation ist der linke Rand fixiert.
Die anderen Ränder sind frei bezüglich der passiven Traktion.

## Aufbau

### Symbolische Beschreibung

- R/expression_tree.R
  Grundstruktur mathematischer Ausdrücke:
  Felder, Operationen, Argumente und Tensorformen.

- R/operators.R
  Symbolische Operatoren zum Aufbau der Gleichungen,
  beispielsweise Gradient, Divergenz und Tensor-Kontraktion.

### Geometrie und physikalische Felder

- R/boundary_conditions.R
  Gitterpunkte, Randklassifikation und fixierte Zustände.

- R/fields.R
  Material-, Dichte- und Kraftfelder, Aktivierung sowie
  Zuordnung numerischer Feldwerte zu symbolischen Koeffizienten.

### Übersetzung in Differentialgleichungen

- R/discretization.R
  Übersetzt symbolische Ausdrücke an Gitterpunkten in
  räumlich diskretisierte Ausdrücke.

- R/ode_system.R
  Baut daraus die Zustandsliste und die rechten Seiten
  des ODE-Systems auf. Fixierte Zustände werden entfernt.

- R/equation_generator.R
  Zentraler Einstieg zur Erzeugung des Kelvin-Voigt-Modells.

### Beispiele, Simulation und Tests

- R/model_script.R
  Beispielkonfiguration zur Modellerzeugung mit einfachen Prüfungen.

- R/dmod2_model.R
  Übergabe an dMod2, Kompilierung und Simulation.

- R/test_simulation.R
  Kontraktionssimulation mit Kraftkalibrierung, Diagnostik,
  PDF-Ausgaben und GIF-Export.

- R/test_active_forces.R
  Tests der aktiven Kraftkonstruktion.

- R/Test_analytic_solution.R
  Vergleich mit einer analytischen Referenzlösung.

## Ablauf

Symbolische Felder und Operatoren
→ Modellgleichung
→ räumliche Diskretisierung und Randbedingungen
→ ODE-System
→ dMod2-Simulation
→ Auswertung und Animation

## Beispiel starten

Aus dem Projektverzeichnis:

    source("R/test_simulation.R")

Benötigt werden dMod2 und für den GIF-Export zusätzlich FFmpeg.

Die aktuelle angepasste Simulation verwendet 40 × 15 Knoten.
Gittergröße, Materialparameter, Dämpfung und Aktivierung werden
am Anfang von test_simulation.R eingestellt.

Die Ausgabe wird im Verzeichnis output/ gespeichert.

## Hinweise

- Alle physikalischen Größen müssen konsistente Einheiten verwenden.
- Beispielparameter sind nicht automatisch biologisch kalibriert.
- Mehr Knoten bei gleichem Abstand vergrößern das Gebiet.
- Ausgabezeitpunkte und interne Integrationsschritte sind verschieden.
- Eine korrekte numerische Referenzlösung ersetzt keine
  Validierung anhand experimenteller Daten.
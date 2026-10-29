# Leseplan: von der Abstraktion zur Simulation

## Ebene 1 – Symbolische Sprache
01. R/expression_tree.R
    Datenstruktur für mathematische Ausdrücke.

02. R/operators.R
    Operationen auf diesen Ausdrücken.

## Ebene 2 – Geometrie und physikalische Felder
03. R/boundary_conditions.R
    Gitter und Randbedingungen.

04. R/fields.R
    Materialfelder, Parameter und aktive Kräfte.

## Ebene 3 – Übersetzung in ein numerisches Modell
05. R/discretization.R
    Räumliche Diskretisierung symbolischer Ausdrücke.

06. R/ode_system.R
    Aufbau der Zustände und Differentialgleichungen.

07. R/equation_generator.R
    Zusammensetzen des Kelvin-Voigt-Modells.

## Ebene 4 – Anwendung
08. R/model_script.R
    Beispiel für die Modellerzeugung.

09. R/dmod2_model.R
    Anbindung an den numerischen Solver.

10. R/test_simulation.R
    Kontraktionssimulation, Auswertung und Animation.

## Begleitende Tests
- R/test_active_forces.R
- R/Test_analytic_solution.R

Hinweis:
Die Nummerierung ist eine didaktische Lesereihenfolge,
keine vorgeschriebene source()-Reihenfolge.
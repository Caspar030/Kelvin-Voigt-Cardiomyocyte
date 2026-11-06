# Offene Fragen zur Materialmodellierung



ToDo's

1. Festplatte nehmen zu Viv

2. Für die Präsentation der Doktorarbeit eine Präsentation machen mit Animationen wie bei 3blue1brown

2. Beliebige Zusammenhängende Mengen erlauben. Nicht nur noch gitter.
3. Boundary conditions einarbeiten: Glass rods. Ein glass rod im inneren: fixierte Punkte. 
Ein anderer Glass rod fixiert nur knoten miteinander fest, lokal.

4. Welche Parameter haben den größten Einfluss darauf, dass sich die Kontraktionskraft ändert
- veränderte Referenzkonfiguration: zueinander verschobene Sarkomere?
- Steifigkeit?

5. Bei Raffi sitzen.





## Soll sich die Steifigkeit während der Kontraktion ändern?

Aktuell wird der Elastizitätstensor E(X) an den Referenzpositionen
ausgewertet und bleibt während der Simulation unverändert.
Die Verformung verändert bereits die elastische Spannung;
eine veränderliche Steifigkeit wäre eine zusätzliche Modellannahme.

| Abhängigkeit | Bedeutung |
| --- | --- |
| E(X) | Unterschiedliche Steifigkeit an verschiedenen Materialorten |
| E(X, ε) | Material wird bei Verformung steifer oder weicher |
| E(X, a(t)) | Steifigkeit verändert sich mit der Aktivierung |
| E(X, t) | Explizite zeitliche Veränderung, beispielsweise durch Umbau |

Dabei bezeichnet X die Referenzposition, ε die Verzerrung
und a(t) den zeitlichen Aktivierungsgrad.

## Zu klären

- Reicht ein zeitlich konstantes E(X) für unsere Fragestellung aus?
- Gibt es experimentelle Hinweise auf verformungs- oder
  aktivierungsabhängige Steifigkeit?
- Ist die lineare Kleinverformungsannahme bei etwa 15 %
  lokaler Verkürzung ausreichend?
- Wie müssten Materialgesetz und Implementierung für eine
  zustandsabhängige Steifigkeit erweitert werden?

Wichtig: E(X + u) beschreibt nicht automatisch eine
verformungsabhängige Steifigkeit. Die Auswertung eines festen
räumlichen Profils an der verschobenen Position ist von einer
Abhängigkeit von der lokalen Verzerrung zu unterscheiden.
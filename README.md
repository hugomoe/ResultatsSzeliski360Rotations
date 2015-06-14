# ResultatsSzeliski360Rotations

Le code pour tester 360 rotations (ou moins) qui retombent sur l'identité.
Le fichier run360.exe, après compilation, prend en argument une image , la fait tourner N fois de 2pi/N radians et donne le résultat de l'expérience (différence L2, L1 et Linfini à l'originale, temps d'exécution). Il crée une image modifiée, en indiquant dans son titre le type d'interpolation, les éventuels paramètres, et si N!=360, le nombre N de rotations. Il stocke les résultats numériques dans un fichier du même nom (remplaçant le précédent fichier du même nom) et à la suite du fichier all_result.c (le créant si besoin).

Pour lancer des tests, renommer d'abord le fichier "all_results.txt" afin que le programme n'écrive pas à la suite de celui existant et qu'on ne mélange pas les résultats obtenus sur différents ordinateurs.

Pour changer de type d'interpolation ou un paramètre d'interpolation (ordre de la B-Spline, beta du raised-cosine), modifier l'en-tête du fichier "run360.c" (entre "///choose an interpolation" et "///end of the choice"). Recompiler ensuite.


Attention, le code tel quel imprime les images en .png. Ceux (celui) qui n'ont pas accès aux bibliothèques libpng et/ou qui ne peuvent pas écrire en .png devront se passer du code, ou le modifier en conséquences.

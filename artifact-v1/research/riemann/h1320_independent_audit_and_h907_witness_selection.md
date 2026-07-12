# H1320 — Audit indépendant et choix du témoin canonique pour H907

Classification: `h1320_exact_math_pass_certificate_float_bug_h1320_canonical_witness`

## Verdict court

Les formules mathématiques de H1320 sont exactes. Une reconstruction autonome
avec `fractions.Fraction`, sans importer le générateur H1320, retrouve :

```text
q2,q3,q4       = 1/2, 1/2, 5/8,
a0,...,a4      = 1, 1, 1/2, 1/8, 5/256,
T3/T2=theta2   = 5/12,
D2             = 13/64,
Disc J^(3,0)   = -27/16,
det Toeplitz   = -5/256.
```

Les identités rationnelles H803/H804/D ont aussi été revérifiées directement
pour `3<=n<=199`; leurs formes symboliques rendent la preuve valable pour tout
`n>=3`.

H1320 doit être retenu comme **témoin minimal canonique** pour l'invalidation du
pont générique H1319/H920 vers H907.

## Vérification des quotients

La construction est

```text
q2 = 1/2,
q_n = (2n-3)/(2n),  n>=3.
```

Pour `n>=3`, calcul direct :

```text
q_{n+1}-q_n = 3/(2n(n+1)) > 0.
```

Comme `q_n<1`, les rapports `R_n=a_n/a_{n-1}` décroissent. Leur produit se
ferme en

```text
R_n = 2 binom(2n-2,n-1)/(n 4^(n-1)).
```

La borne centrale `binom(2m,m)<=4^m` donne `R_n<=2/n`, puis

```text
a_n <= 2^(n-1)/n!.
```

La série `sum a_n z^n` est donc entière : le témoin n'est pas seulement un
préfixe fini artificiel.

## Vérification H803 et H804

Avec `C_n=(2n)(2n-1)` et

```text
theta_n=C_{n-1}C_{n+1}/C_n^2,
```

l'identité entre les quotients de `a_n` et ceux des moments
`M_n=(2n)!a_n` est

```text
T_{n+1}/T_n = theta_n q_{n+1}/q_n.
```

Pour H1320, `n>=3` :

```text
T_{n+1}/T_n
  = (n-1)(2n+1)/(n(2n-1)),

T_{n+1}/T_n - theta_n
  = 3(n-1)(2n+1)/(n^2(2n-1)^2) > 0,

T_{n+1}/T_n - (1-1/n^2)
  = (n-1)/(n^2(2n-1)) > 0.
```

Au cas de base `n=2`, la reconstruction

```text
M0,M1,M2,M3 = 1,2,12,90
```

donne exactement

```text
T3/T2 = 5/12 = theta2.
```

Ainsi H803 tient pour tout `n>=2`, et la barrière unitaire H804 tient pour tout
`n>=3`. Le fait qu'elle échoue à `n=2` est sans incidence : H804 est une
barrière éventuelle, les petits indices relevant de H803 exact.

## Vérification des gardes D

En posant `e_n=1-q_n`, le garde utilisé dans H851/H877 est

```text
D_n=e_{n+1}^2-q_{n+1}^2 e_n e_{n+2}.
```

Le calcul donne

```text
D2 = 13/64 > 0,

D_n = 9(12n-1)/(16n(n+1)^2(n+2)) > 0,  n>=3.
```

Ajouter les gardes D translatés au paquet H803/H804 ne supprime donc pas le
contre-modèle.

## Deux falsificateurs indépendants

Avec `gamma_j=j!a_j`, les quatre premiers coefficients sont

```text
gamma = (1,1,1,3/4).
```

Donc

```text
J_gamma^(3,0)(x)=1+3x+3x^2+(3/4)x^3
```

a pour discriminant exact

```text
-27/16<0.
```

Il possède une paire de racines non réelles.

Indépendamment, le mineur Toeplitz

```text
rows=(0,1,2,3), cols=(1,2,3,4)
```

vaut exactement

```text
-5/256<0.
```

Le premier falsificateur bloque directement l'hyperbolicité Jensen ; le second
bloque PF4 et donc PF-infinity sans utiliser de calcul de racines.

## Comparaison avec les deux témoins H907 supplémentaires

| Témoin | Antécédents réellement certifiés | Échec | Rôle recommandé |
| --- | --- | --- | --- |
| H1320 | suite entière infinie, H803 tout `n>=2`, H804 tout `n>=3`, `q` croissant, PF2, tous les gardes D | Jensen degré 3 et PF4 | **canonique** pour réfuter le pont H920 vers H907 |
| `h907_pf3_to_pf4_exact_countermodel` | fenêtre rationnelle `0..8`, PF3 complet, `q` croissant | PF4 | témoin auxiliaire de la stricte hiérarchie PF |
| `h907_jensen_degree4_exact_countermodel` | fenêtre `0..8`, PF4, tous Jensen degrés 2–3, `q` croissant | Jensen degré 4 | témoin complémentaire si l'on ajoute les petits degrés comme hypothèses |

H1320 est le plus fort pour la question réellement posée parce qu'il satisfait
le paquet exact et global fermé par H1319, avec une formule valable à l'infini.
Il est aussi le plus petit algébriquement : quatre coefficients suffisent au
falsificateur Jensen et cinq au falsificateur Toeplitz.

Le témoin Jensen degré 4 reste utile contre une réparation particulière :
« ajoutons les degrés 2 et 3 comme base d'induction ». Il montre que cette base
finie ne fournit toujours aucune induction abstraite. Mais sa nature de fenêtre
finie le rend moins canonique que H1320.

## Caveat d'implémentation du certificat corrigé

Deux validations passent :

```text
python tools/rh_h1320_factorial_barrier_all_degree_counterexample_canonical.py --check-report
python tools/rh_h1320_factorial_barrier_independent_audit.py
```

En revanche, au moment de cet audit :

```text
python tools/rh_h1320_factorial_barrier_all_degree_certificate.py --check-report
```

échoue. La cause est locale au code : `theta2` est calculé avec la division
Python `/` avant conversion SymPy, ce qui produit

```text
0.416666666666667
```

au lieu de l'exact `5/12`. Par conséquent le calcul frais marque
`n2_factorial_ratio_reconstructed=False`, tandis que le JSON épinglé contient
la valeur rationnelle correcte.

La réparation doit remplacer cette division par `sp.Rational(num,den)`. Ce bug
n'invalide pas H1320 : la reconstruction Fraction autonome, l'ancien certificat
canonique et toutes les identités symboliques confirment `theta2=5/12`.

## Recommandation H907

Utiliser H1320 comme unique contre-modèle principal et formuler désormais le
verrou ainsi :

```text
Le paquet Xi H803/H804/H920, même renforcé par PF2 et les gardes D,
ne contient aucun mécanisme tous degrés.

Toute fermeture H907 doit introduire un invariant Xi-spécifique uniforme en d
qui exclut explicitement la suite H1320.
```

Conserver le témoin Jensen degré 4 en annexe comme test anti-réparation pour
toute proposition qui n'ajouterait que les degrés 2 et 3. Ne pas poursuivre de
nouveaux scans PF4/PF5 avant qu'un invariant Xi-spécifique précis ne soit
énoncé.


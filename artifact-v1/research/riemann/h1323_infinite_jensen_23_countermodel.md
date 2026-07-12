# H1323 — Contre-modèle infini exact au patch « Jensen degrés 2–3 »

Classification: `h1323_infinite_exact_countermodel_to_low_degree_jensen_patch`

## Verdict

Le contre-modèle fini de `h907_jensen_degree4_exact_countermodel.md` admet un
prolongement rationnel positif à tous les indices. La suite prolongée vérifie
simultanément :

```text
- PF2 sur toute la matrice Toeplitz unilatérale ;
- q_(n+1) >= q_n et q_n<1 pour tout n>=2, donc H803 global ;
- q_2=13/32<1/2 ;
- la barrière unitaire H804 pour tout n>=8 ;
- D_n>0 pour tout n>=2, où D_n est le déficit traduit H851 ;
- Phi(z)=sum a_n z^n entière ;
- tous les polynômes de Jensen de degrés 2 et 3, à tout décalage, strictement
  hyperboliques ;
```

mais `J_gamma^(4,0)` possède une paire de racines non réelles.

Cela invalide donc un patch sensiblement plus fort que « ajoutons les degrés
2–3 » : même toutes les hypothèses scalaires actuellement disponibles dans la
route H811/H851, plus H804 dans la queue, ne permettent pas de passer au degré
4 par un argument universel.

Ce résultat ne porte pas sur les vrais coefficients de Xi. Il montre que le
pont manquant doit contenir une information propre à Xi et uniforme en degré.

## 1. Construction infinie

Prendre

```text
(q_2,...,q_8)
 = (13/32, 141/256, 85/128, 89/128,
    193/256, 201/256, 201/256)
```

et, pour tout `n>=8`, poser

```text
c   = 201/208,
q_n = c (2n-3)/(2n).
```

Le raccord est exact car

```text
c (2*8-3)/(2*8) = (201/208)(13/16) = 201/256 = q_8.
```

Définir ensuite

```text
R_1=a_0=a_1=1,
R_n/R_(n-1)=q_n,
a_n/a_(n-1)=R_n,
gamma_n=n! a_n.
```

Tous ces nombres sont rationnels et strictement positifs.

## 2. PF2, H803 et H804

Les différences dans le préfixe sont

```text
q_3-q_2, ..., q_8-q_7
 = 37/256, 29/256, 1/32, 15/256, 1/32, 0.
```

Dans la queue,

```text
q_(n+1)-q_n = 603/[416 n(n+1)] > 0,
0<q_n<c=201/208<1.
```

Ainsi `q_n` est non décroissante sur tout `n>=2`, ce qui donne H803. De plus,

```text
R_(n+1)/R_n=q_(n+1)<1,
```

donc les rapports `R_n=a_n/a_(n-1)` décroissent. La suite positive `a_n` est
log-concave sans trou interne, ce qui est exactement PF2 pour la matrice
Toeplitz unilatérale.

Pour H804, avec

```text
C_n=(2n)(2n-1),
tau_n=q_n C_n/C_(n-1),
```

le calcul exact dans la queue donne, pour tout `n>=8`,

```text
tau_(n+1)/tau_n
 = (n-1)(2n+1)/[n(2n-1)],

tau_(n+1)/tau_n - (1-1/n^2)
 = (n-1)/[n^2(2n-1)] > 0.
```

La barrière unitaire est donc satisfaite strictement.

## 3. La série génératrice est entière

Comme `q_n<c<1` pour tout `n>=2`,

```text
R_n <= c^(n-1),
a_n <= c^[n(n-1)/2],
a_n^(1/n) <= c^[(n-1)/2] -> 0.
```

Le critère de la racine montre que

```text
Phi(z)=sum_(n>=0) a_n z^n
```

est entière.

## 4. Le déficit traduit H851 est globalement positif

Poser

```text
D_n=(1-q_(n+1))^2-q_(n+1)^2(1-q_n)(1-q_(n+2)).
```

Les six valeurs de raccord sont

```text
D_2 = 37926823/268435456,
D_3 = 28183907/536870912,
D_4 = 28382139/536870912,
D_5 = 50157087/2147483648,
D_6 = 58256935/4294967296,
D_7 = 142163945/6979321856.
```

Elles sont toutes strictement positives. Pour `n>=8`, la substitution de la
formule de queue se factorise en

```text
D_n =
 (2244592 n^4
  + 297357088 n^3
  + 13154707832 n^2
  + 192552254264 n
  - 15372297693)
 /
 (29948379136 n (n+1)^2 (n+2)).
```

Le dénominateur est positif. Dans le numérateur, tous les coefficients non
constants sont positifs et déjà

```text
192552254264 n - 15372297693 > 0
```

pour `n>=1`. Par conséquent `D_n>0` pour tout `n>=2`.

### Portée exacte vis-à-vis de H811

Le témoin satisfait donc les hypothèses scalaires utilisées par H811 :

```text
positivité + PF2,
q_2<=1/2,
q croissant,
déficits traduits D_n>=0.
```

On peut donc lui appliquer la conclusion démontrée par H811 : les mineurs PF3
à lignes consécutives sont non négatifs. Cela ne prouve pas PF3 complet : les
mineurs à lignes clairsemées restent hors de la portée de H811, exactement
comme l'indique son énoncé. Aucune revendication PF3 globale n'est faite ici.

## 5. Jensen degrés 2 et 3 à tous les décalages

Écrire

```text
P_m=gamma_m/gamma_(m-1),
Q_m=P_m/P_(m-1)=m q_m/(m-1).
```

Après division par `gamma_n` et le changement positif
`y=P_(n+1)X`, les polynômes de Jensen normalisés deviennent

```text
d=2 : 1 + 2y + u y^2,
d=3 : 1 + 3y + 3u y^2 + u^2 v y^3,

u=Q_(n+2), v=Q_(n+3).
```

Pour les décalages `n=0,...,5`, les discriminants normalisés de degré 2 sont

```text
3/4, 89/128, 11/24, 67/128, 61/160, 43/128,
```

et ceux de degré 3 sont

```text
6829930341/17179869184,
9833091326901/70368744177664,
119928981575/824633720832,
96969867412005/1125899906842624,
2602796310705933/43980465111040000,
4456654550901/70368744177664.
```

Ils sont strictement positifs.

Pour tout décalage `n>=6`, les deux quotients nécessaires sont déjà dans la
queue :

```text
u = c(2n+1)/[2(n+1)],
v = c(2n+3)/[2(n+2)].
```

Le discriminant quadratique normalisé vaut

```text
4(1-u) = (14n+215)/[104(n+1)] > 0.
```

Pour le cubique, poser

```text
F(u,v)=u^2 v^2-6uv+4u+4v-3.
```

Son discriminant normalisé est `-27u^2 F`. Or le calcul exact donne

```text
F = -3
    (1509200 n^4
     + 73874752 n^3
     + 1236652088 n^2
     + 7568113136 n
     + 5859746333)
    /
    (29948379136 (n+1)^2 (n+2)^2)
  < 0.
```

Le discriminant cubique est donc strictement positif. Tous les Jensen de
degrés 2 et 3 sont ainsi hyperboliques pour tous les décalages, sans troncature
numérique.

## 6. Échec exact au degré 4

Les quatre premiers quotients ne sont pas modifiés par le prolongement. On
retrouve

```text
J_gamma^(4,0)(X)
 = 1
   + 4 X
   + (39/8) X^2
   + (71487/32768) X^3
   + (11138032035/34359738368) X^4.
```

Son discriminant exact est

```text
-691976798016209978769928197
/158456325028528675187087900672
< 0.
```

Le quartique réel possède donc exactement une paire conjuguée non réelle. Le
degré 4 n'est pas hyperbolique.

## 7. Conclusion logique

Le contre-modèle invalide toute implication universelle de la forme

```text
PF2
+ q croissant et q_2<=1/2
+ H804 dans la queue
+ tous les D_n traduits positifs
+ génératrice entière positive
+ Jensen hyperbolique aux degrés 2 et 3 pour tout décalage
=> Jensen hyperbolique au degré 4.
```

Il invalide en particulier l'idée de combler H907 en ajoutant simplement les
degrés 2 et 3 à H803/H851. Le prochain lemme utile doit être véritablement
tous degrés ou doit introduire une identité intégrale/déterminantale propre au
noyau Xi qui exclut ce témoin.

## 8. Artefacts stricts

```text
tools/rh_h1323_infinite_jensen_23_countermodel_canonical_certificate.py
research/riemann/h1323_infinite_jensen_23_countermodel.json
tools/rh_h1323_infinite_jensen_23_countermodel_independent_audit.py
research/riemann/h1323_infinite_jensen_23_countermodel_independent_audit.json
```

Les commandes strictes sont

```text
python tools/rh_h1323_infinite_jensen_23_countermodel_canonical_certificate.py --check-report
python tools/rh_h1323_infinite_jensen_23_countermodel_independent_audit.py --check-report
```

Elles terminent toutes deux avec le code `0` et `all_checks_pass=true`. Le
second vérificateur n'importe ni SymPy ni le code primaire : il reconstruit les
identités rationnelles avec `fractions.Fraction` et une arithmétique
polynomiale univariée indépendante.

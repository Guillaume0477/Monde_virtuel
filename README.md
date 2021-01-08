# Monde_virtuel
# TP 2 relatif à la gestion de carte de hauteurs et création de caractéristiques
# Updated : TP3 : ajout de la végétation dont la distribution dépend des caractéristiques de terrain
# Auteurs : BRALET Antoine & DURET Guillaume

Le code permettant de créer les différentes caractéristiques des cartes de hauteurs 
et d'enregistrer différentes cartes (slope, slopeAVG, stream area, stream power, 
accessibility,... se trouve dans le dossier TP2. 

Updated : Ce même dossier contient également les classes et fonctions permettant 
de créer une distribution de végétaux à partir de l'humidité, le stream power, 
l'accessibilité et la pente du terrain. Plusieurs images sont fournis dans le 
dossier /ImagesToTest afin de tester ce code sur différentes tailles de terrain. 

Le code permettant de visualiser la carte de hauteur en 3D est disponible dans le
dossier tp_parallax_mapping. Il suffit de modifier les lignes d'appels aux 
textures pour pouvoir visualiser une carte donnée sur le terrain généré par la carte
de hauteur ou le fichier .obj créé précédemment. Ce code a été implémenté lors d'un 
précédent TP à CPE Lyon dans le cours de T. DUPONT. Nous l'avons donc remanié pour 
qu'il soit utilisable ici.

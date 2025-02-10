# Projekt zawiera:

**Shadow Mapping** – wykorzystująca mapę głębokości, aby określić, które fragmenty sceny są zacienione względem źródła światła.

**Tekstury, normal mapping** - obrazy nakładane na modele 3D w celu nadania im realistycznego wyglądu, np. kolorów, nierówności, wzorów czy szczegółów powierzchni. Normal mapping polepsza też realizm oświetlenia.

**Algorytm Boids** – symulacja ruchu i zachowania stada ptaków, algorytm uwzględnia zasady separacji, dopasowania i spójności.

**Detekcja kolizji AABB (Axis-Aligned Bounding Box)** – pozwala sprawdzić, czy obiekty nachodzą na siebie na podstawie zdefiniowanych ograniczających je „pudełek”, co umożliwia boidom (ptakom) unikanie przeszkód.

Boidsy omijają gdy, mają kolizję z pokojem, drzewami. Nie wlatują również w teren.

**Generowanie Terenu** – zadanie proceduralne, przy użyciu algorytmu szumu wraz ze smoothingiem i interpolacją. Wysokość generowana jest za pomocą generateHeight(), która wykorzystuje dodatkowo oktawy i roughness w celu zwiększenia “randomowości” terenu. Na bazie tej funkcji stawiane są później drzewa (.obj), które dzięki wynikowi znajdują się w naturalnie wyglądających miejscach.

**Skybox** – oteksturowana Cubemapa wraz z dedykowanymi shaderami skybox.vert i skybox.frag imitująca odległe tereny

# Sterowanie projektem:

**Przyciski C i V** włączają i wyłączają shadow mapping.

**Przyciski P i O** zatrzymują i ruszają ptaki.

**Klawisze 1 i 2** przyciemniają i rozjaśniają scenę.

**Klawisze K i L** ustawiają poziom latania boidów. W górę lub w dół.

# Podział pracy:

**Paulina Kaczmarek** : algorytm boids, interakcja – zatrzymywanie boidów, ustawianie poziomu wysokości latania boidsów, detekcja kolizji AABB

**Paulina Śmierzchalska** : tekstury, normal mapping, shadow mapping, interakcja – wyłączanie i włączanie shadow mappingu

**Michał Hołyniewski** : proceduralne generowanie terenu, skybox, wykończenie wizualne

# Wyzwania:

- Cieniowanie na terenie, samo shadow mapping
- Dobranie odpowiednich parametrów do budowy bouding box
- Naprawa rozbieżności funkcji wysokości terenu a obiektów i ujednolicenie systemu współrzędnych

Jeśli pojawia się błąd z bibliotekami, wystarczy podmienić biblioteki z ćw9.



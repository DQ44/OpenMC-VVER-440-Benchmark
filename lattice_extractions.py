"""
Visualisation d'un openmc.HexLattice a partir de lattice.indices et
lattice.universes uniquement (pas de fonction de plot OpenMC).

Orientation de l'hexagone affiche (verifie visuellement, pas de trou ni
chevauchement) :
    - core.orientation == 'x'  -> hexagone pointe en haut / en bas
    - core.orientation == 'y'  -> hexagone plat en haut / en bas (defaut)
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import RegularPolygon
from matplotlib.lines import Line2D
from math import sqrt


def hex_lattice_positions(lattice):
    """
    Calcule la position cartésienne (X, Y) du centre de chaque élément du
    lattice, indexée par le même (r, i) que lattice.indices et
    lattice.universes[r][i].

    Retour
    ------
    dict : {(r, i): (X, Y)}
    """
    num_rings = lattice.num_rings
    orientation = lattice.orientation          # 'x' ou 'y'
    pitch = lattice.pitch[0]
    cx, cy = lattice.center[0], lattice.center[1]

    def get_universe_index(x, a):
        # Reproduction de HexLattice.get_universe_index (openmc/lattice.py)
        z = -a - x
        g = max(abs(x), abs(a), abs(z))
        i_ring = num_rings - 1 - g
        if x >= 0:
            i_within = x if a >= 0 else 2 * g + z
        else:
            i_within = 3 * g - x if a <= 0 else 5 * g - z
        if orientation == 'x' and g > 0:
            i_within = (i_within + 5 * g) % (6 * g)
        return i_ring, i_within

    positions = {}
    for a in range(-num_rings + 1, num_rings):
        for x in range(-num_rings + 1, num_rings):
            z = -a - x
            g = max(abs(x), abs(a), abs(z))
            if g >= num_rings:
                continue  # indice hors du lattice (is_valid_index -> False)

            r, i = get_universe_index(x, a)

            # Reproduction de HexLattice.get_local_coordinates, en résolvant
            # pour le centre de l'élément (point = centre => coord locale 0)
            if orientation == 'x':
                X = cx + (x + 0.5 * a) * pitch
                Y = cy + sqrt(0.75) * pitch * a
            else:  # orientation 'y', valeur par défaut d'OpenMC
                X = cx + sqrt(0.75) * pitch * x
                Y = cy + (0.5 * x + a) * pitch

            positions[(r, i)] = (X, Y)

    return positions


def plot_hex_core(lattice, colors, ax=None, figsize=(11, 11)):
    """
    Trace le coeur hexagonal, sans axes, avec une légende a hexagones.

    Paramètres
    ----------
    lattice : openmc.HexLattice
    colors : dict {universe.id: (couleur_matplotlib, label)}
        Une entrée par type d'universe (assembly16, assembly24, assembly425,
        moderator_universe, ...). C'est cette couleur qui est utilisée à la
        fois pour les hexagones et pour la légende.
    """
    positions = hex_lattice_positions(lattice)
    pitch = lattice.pitch[0]
    radius = pitch / sqrt(3)  # rayon centre -> sommet d'un hexagone régulier

    # 'x' -> hexagone pointe en haut/bas ; 'y' -> hexagone plat en haut/bas
    mpl_orientation = 0.0 if lattice.orientation == 'x' else np.pi / 6

    if ax is None:
        _, ax = plt.subplots(figsize=figsize)

    for (r, i) in lattice.indices:
        univ = lattice.universes[r][i]
        X, Y = positions[(r, i)]
        color, _ = colors[univ.id]
        ax.add_patch(RegularPolygon(
            (X, Y), numVertices=6, radius=radius,
            orientation=mpl_orientation,
            facecolor=color, edgecolor='black', linewidth=0.5
        ))

    ax.set_aspect('equal')
    ax.autoscale_view()
    ax.axis('off')  # coeur propre, sans axes ni graduations

    # Légende : petits hexagones a bord noir fin, meme couleur que le graphique
    handles = [
        Line2D([0], [0], marker='H', linestyle='', markersize=16,
               markerfacecolor=c, markeredgecolor='black', markeredgewidth=0.8,
               label=label)
        for (c, label) in colors.values()
    ]
    ax.legend(handles=handles, loc='center left', bbox_to_anchor=(1.0, 0.5),
              frameon=False)

    return ax


# ---------------------------------------------------------------------
# Exemple d'utilisation avec votre géométrie (à adapter / décommenter)
# ---------------------------------------------------------------------
if __name__ == '__main__':
    # import openmc
    # ... (définition de moderator_universe, assembly16, assembly24,
    #      assembly425, core comme dans votre script) ...

    colors = {
        # moderator_universe.id: ('#cfe8ff', 'Moderateur'),
        # assembly16.id:         ('#4caf50', 'Assembly 1.6%'),
        # assembly24.id:         ('#ff9800', 'Assembly 2.4%'),
        # assembly425.id:        ('#e53935', 'Assembly 4.25%'),
    }

    # fig, ax = plt.subplots(figsize=(12, 12))
    # plot_hex_core(core, colors, ax=ax)
    # plt.tight_layout()
    # plt.savefig('core_layout.png', dpi=150, bbox_inches='tight')
    # plt.show()
    pass

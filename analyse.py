#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  analyse.py
#
#  Copyright 2019 Paulo Mappa <paulo@ile>

"""
------------------------------------------------
Unidade Varginha do CEFET-MG
Departamento de Computação e Engenharia Civil
Graduação em Engenharia Civil
Disciplina Análise Matricial de Estruturas
Prof. Paulo Mappa
2019/2
-------------------------------------------------

Análise linear de estruturas reticuladas.

"""


def desenha(inci, coord, titulo='     ', cor='r'):
    """
    Desenha o reticulado a a partir das coordenadas nodais(coord) e
    incidências dos lementos(inci)
    titulo: string com o título que aparecerá no desenho da estrutura
             por padrão nehum título ararecerá.
    cor: Cor usdad nas linhas do reticulado. Padrão Matiplotlib, ou seja:
        g: green, b:black, r:red ...
    """
    import matplotlib.pyplot as plt
    elemento = 0
    for i in inci:
        elemento = elemento + 1
        no_inicial = i[0]
        no_final = i[1]
        xi = coord[no_inicial-1, 0]
        xf = coord[no_final-1, 0]
        yi = coord[no_inicial-1, 1]
        yf = coord[no_final-1, 1]
        plt.plot((xi, xf), (yi, yf), cor)
        plt.text(xi, yi, no_inicial)
        plt.text(xf, yf, no_final)
        ele = "[{}]".format(elemento)  # a chave é o coringa!
        plt.text(0.90*xi+(xf-xi)/2, yi+(yf-yi)/2, ele)
    plt.title(titulo)
    plt.show()


def dados(problema=1):
    """
    Entradas de dados.
        deve ser editada com os dados da estrutura:
            coord: coordenadas nodais;
            inci: incidências dos elementos;
            rest: restrições nodais;
            e: módulo de elasticidade;
            a: área dos elementos.
            f: açoes externas nos nós
        (e, a, coord, inci, rest, f) = dados(problema)
            dados                 descrição                referência
              1        trelica plana de três barras          pag 49, notas
                                                            Prof Ronaldo Batista
    """
    import numpy as np
    e = 200e6
    a = 0.0015
    if problema == 1:
        coord = np.array([[0, 0], [2, 0], [0, 1.5]], float)
        inci = np.array([[1, 2], [2, 3], [1, 3]], int)
        rest = np.array([[0, 1], [1, 1], [0, 0]], int) # 0: livre;
                                                       # 1: restringido
        f = np.array([[0.0, -20.0, 0.0, -45.0, 0.0, 0.0]], float)
    return (e, a, coord, inci, rest, f)


def trelica(e, a, no_i, no_j, l, s, c):
    """
    Monta a matriz do elemnto de treliça plana nop referencial global.
        no_i: nó inicial do elemento;
        no_j: nó final do elemento;
        s,c: seno e cosseno diretor do elemento;
        l: comprimento do elemento;
        k: matriz de rigidez do elemento - referêncial Global;
        endereco_global: contém as direções globais correspondentes às direçoes
                         das deslocabilidades locais do elemento.
        (k, endereco_global) = trelica(e, a, no_i, no_j, l, s, c)
    """
    import numpy as np
    k = np.zeros((4, 4))
    k[0, 0] = e*a/l
    k[0, 2] = -e*a/l
    k[2, 2] = k[0, 0]
    k[2, 0] = k[0, 2]
    r = np.array([[c, s, 0, 0], [-s, c, 0, 0], [0, 0, c, s], [0, 0, -s, c]],
                 float)
    k = np.dot(k, r)
    k = np.dot(r.transpose(), k)
    endereco_global = np.array([[2*no_i - 2, 2*no_i - 1, 2*no_j - 2, 2*no_j - 1
                                 ]], int)
    return(k, endereco_global)


def monta_matriz_global_trelica(k_elemento, endereco_global, k_estrutura):
    """
    Monta a matriz de rigidez da treliça no referencial global.
    k_elemento: matriz de rigidez do elemento no referencial global.
    k_estrutura: Matriz de rigidez da treliça no referencial global.

    """
    import numpy as np
    #  import pdb        # depurador
    #  pdb.set_trace()   # depurador para aqui
    k_estrutura[endereco_global[0, 0], endereco_global[0, 0]] =\
        k_estrutura[endereco_global[0, 0], endereco_global[0, 0]] +\
        k_elemento[0, 0]

    k_estrutura[endereco_global[0, 0], endereco_global[0, 1]] =\
        k_estrutura[endereco_global[0, 0], endereco_global[0, 1]] +\
        k_elemento[0, 1]

    k_estrutura[endereco_global[0, 1], endereco_global[0, 0]] =\
        k_estrutura[endereco_global[0, 0], endereco_global[0, 1]]

    k_estrutura[endereco_global[0, 0], endereco_global[0, 2]] =\
        k_estrutura[endereco_global[0, 0], endereco_global[0, 2]] +\
        k_elemento[0, 2]

    k_estrutura[endereco_global[0, 2], endereco_global[0, 0]] =\
        k_estrutura[endereco_global[0, 0], endereco_global[0, 2]]

    k_estrutura[endereco_global[0, 0], endereco_global[0, 3]] =\
        k_estrutura[endereco_global[0, 0], endereco_global[0, 3]] +\
        k_elemento[0, 3]

    k_estrutura[endereco_global[0, 3], endereco_global[0, 0]] =\
        k_estrutura[endereco_global[0, 0], endereco_global[0, 3]]

    # ---------

    k_estrutura[endereco_global[0, 1], endereco_global[0, 1]] =\
        k_estrutura[endereco_global[0, 1], endereco_global[0, 1]] +\
        k_elemento[1, 1]

    k_estrutura[endereco_global[0, 1], endereco_global[0, 2]] =\
        k_estrutura[endereco_global[0, 1], endereco_global[0, 2]] +\
        k_elemento[1, 2]

    k_estrutura[endereco_global[0, 2], endereco_global[0, 1]] =\
        k_estrutura[endereco_global[0, 1], endereco_global[0, 2]]

    k_estrutura[endereco_global[0, 1], endereco_global[0, 3]] =\
        k_estrutura[endereco_global[0, 1], endereco_global[0, 3]] +\
        k_elemento[1, 3]

    k_estrutura[endereco_global[0, 3], endereco_global[0, 1]] =\
        k_estrutura[endereco_global[0, 1], endereco_global[0, 3]]

    # ---------

    k_estrutura[endereco_global[0, 2], endereco_global[0, 2]] =\
        k_estrutura[endereco_global[0, 2], endereco_global[0, 2]] +\
        k_elemento[2, 2]

    k_estrutura[endereco_global[0, 2], endereco_global[0, 3]] =\
        k_estrutura[endereco_global[0, 2], endereco_global[0, 3]] +\
        k_elemento[2, 3]

    k_estrutura[endereco_global[0, 3], endereco_global[0, 2]] =\
        k_estrutura[endereco_global[0, 2], endereco_global[0, 3]]

    # ---------

    k_estrutura[endereco_global[0, 3], endereco_global[0, 3]] =\
        k_estrutura[endereco_global[0, 3], endereco_global[0, 3]] +\
        k_elemento[3, 3]

    return(k_estrutura)


def reordena(k_estrutura, f_nodais, rest):
    """
    Impoe as condiçoes dos apoios à matriz de rigidez e ao vetor de cargas.

    """

    import numpy as np
    import pdb        # depurador
    correspondencia_aux = np.zeros((2, 2*len(rest)), int)
    correspondencia_aux[0, :] = rest.reshape(1, 2*len(rest))
    for i in range(np.size(correspondencia_aux, 1)):
        if i == 0:
            correspondencia_aux[1, i] = correspondencia_aux[0, i]
        else:
            if correspondencia_aux[0, i] == 0:
                correspondencia_aux[1, i] = 0
            else:
                correspondencia_aux[1, i] = np.max(correspondencia_aux[1,:]) + 1

    correspondencia = np.zeros(np.max(correspondencia_aux), int)
    f = np.zeros(len(correspondencia), float)
    #pdb.set_trace()   # depurador para aqui
    j = 0
    for i in correspondencia_aux[1, :]:
        if i != 0:
            correspondencia[j] = j+1
            j = j + 1
    k = np.zeros((len(correspondencia), len(correspondencia)), int)
    #pdb.set_trace()   # depurador para aqui
    for i in range(len(correspondencia)):
        f[i] = f_nodais[0, correspondencia[i]]
        for j in range(len(correspondencia)):
            k[i, j] = k_estrutura[correspondencia[i], correspondencia[j]]
    #pdb.set_trace()
    return(k, f, correspondencia_aux[1, :])

def deformada(coord, inci, u_restricao, indices):
    """
    Distribui u_resticao, vetor com os deslocamentos nodais, considerando as
    direções generalizadas contidas em indices.  Armazena em indices

    """
    from numpy import reshape, size, zeros
    from analyse import desenha
    import pdb
    aux = reshape(coord, size(coord), 0)
    aux1 = zeros(size(aux))
    j = 0
    for i in indices:
        if i != 0:
            aux1[i] = u_restricao[j]*100
            j = j + 1
    aux = aux + aux1
    k = 0
    for i in range(size(coord,0)):
        for j in range(size(coord,1)):
            coord[i,j] = aux[k]
            k = k + 1
    desenha(inci, coord, titulo='     ', cor='b')


def analyse(problema=1):
    """
    resolve o problema cálculo de deslocamento de uma treliça plana com
    comportamento linear geométrico e físico.

    """
    import pdb        # depurador
    import numpy as np
    (e, a, coord, inci, rest, f_nodais) = dados(problema)
    elemento = 0
    ngl = 2*len(coord)
    k_estrutura = np.zeros((ngl, ngl), float)
    for i in inci:
        elemento = elemento + 1
        no_inicial = i[0]
        no_final = i[1]
        xi = coord[no_inicial-1, 0]
        xf = coord[no_final-1, 0]
        yi = coord[no_inicial-1, 1]
        yf = coord[no_final-1, 1]
        comprimento = ((xf - xi)**2 + (yf - yi)**2)**0.5
        s = (yf - yi)/comprimento
        c = (xf - xi)/comprimento
        (k_elemento_global, endereco_global) =\
            trelica(e, a, no_inicial, no_final, comprimento, s, c)
        k_estrutura = monta_matriz_global_trelica(k_elemento_global,
                                                  endereco_global, k_estrutura)
    (k, f, indices) = reordena(k_estrutura, f_nodais, rest)
    #pdb.set_trace()
    #f = np.reshape(f, size(f), 1)
    u = np.linalg.solve(k, f)
    deformada(coord, inci, u, indices)
    return(u)

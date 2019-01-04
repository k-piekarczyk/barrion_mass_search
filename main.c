//
// Created by Krzysiek on 04.01.2019.
//

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define BUFOR_ROZNICY 99999999

/** Definicja listy wierszy */
typedef struct _wiersz_node {
    double m_barionowa;
    double m_grawitacyjna;
    double m_pedu;
    double R;               // promien chuja
    double omega;           // predkosc katowa
    double P;               // okres rotacji
    double Ro_B;            // centralna gestosc energii
    double N_B;             // gestosc barionowa
    double entalpia;
    double r_eq;            // promien przy rowniku
    double s_r;             // stosunek promieni
    double s_e;             // stosunek energii

    struct _wiersz_node *next;
} WierszNode;

void
appendWiersz(WierszNode **head_ref, double m_barionowa, double m_grawitacyjna, double m_pedu, double R, double omega,
             double P, double Ro_B, double N_B, double entalpia, double r_eq, double s_r, double s_e) {
    WierszNode *new_wiersz = (WierszNode *) malloc(sizeof(WierszNode));
    WierszNode *current_node = *head_ref;

    new_wiersz->next = NULL;
    new_wiersz->m_barionowa = m_barionowa;
    new_wiersz->m_grawitacyjna = m_grawitacyjna;
    new_wiersz->m_pedu = m_pedu;
    new_wiersz->R = R;
    new_wiersz->omega = omega;
    new_wiersz->P = P;
    new_wiersz->Ro_B = Ro_B;
    new_wiersz->N_B = N_B;
    new_wiersz->entalpia = entalpia;
    new_wiersz->r_eq = r_eq;
    new_wiersz->s_r = s_r;
    new_wiersz->s_e = s_e;


    if (*head_ref == NULL) {
        *head_ref = new_wiersz;
        return;
    }

    while (current_node->next != NULL)
        current_node = current_node->next;

    current_node->next = new_wiersz;
}

WierszNode *findByMBar(WierszNode **head_ref, double m_bar) {
    WierszNode *current_node = *head_ref;

    while (current_node != NULL) {
        if (current_node->m_barionowa == m_bar) return current_node;
        current_node = current_node->next;
    }

    return NULL;
}

void destroyList(WierszNode **head_ref) {
    WierszNode *current_node = *head_ref;
    WierszNode *next;

    while (current_node != NULL) {
        next = current_node->next;
        free(current_node);
        current_node = next;
    }

    *head_ref = NULL;
}

/** Przeszukiwanie najbliższych */
typedef struct _para_wierszy {
    WierszNode *mniejszy;
    WierszNode *wiekszy;
} ParaWierszy;

ParaWierszy *znajdzNajblizszaPare(WierszNode **head_ref, double m_bar_szukana) {
    WierszNode *obecny_wezel = *head_ref;

    WierszNode *mniejsza = NULL;
    WierszNode *wieksza = NULL;
    double min_diff = BUFOR_ROZNICY;
    double max_diff = BUFOR_ROZNICY;
    double diff;


    while (obecny_wezel != NULL) {
        diff = fabs(obecny_wezel->m_barionowa - m_bar_szukana);

        if ((obecny_wezel->m_barionowa < m_bar_szukana) && (diff < min_diff)) { // Szukanie mniejszej
            min_diff = diff;
            mniejsza = obecny_wezel;
        } else if ((obecny_wezel->m_barionowa > m_bar_szukana) && (diff < max_diff)) { // Szukanie wiekszej
            max_diff = diff;
            wieksza = obecny_wezel;
        }

        obecny_wezel = obecny_wezel->next;
    }

    if (mniejsza == NULL) {
        printf("Nie znaleziono wiersza z mniejsza masa barionowa niz: %f\n", m_bar_szukana);
        exit(EXIT_FAILURE);
    } else if (wieksza == NULL) {
        printf("Nie znaleziono wiersza z wieksza masa barionowa niz: %f\n", m_bar_szukana);
        exit(EXIT_FAILURE);
    }

    ParaWierszy *paraWierszy = malloc(sizeof(ParaWierszy));
    paraWierszy->mniejszy = mniejsza;
    paraWierszy->wiekszy = wieksza;

    return paraWierszy;
}


int main(int argc, char **argv) {
    int m_barionowa, m_grawitacja;
//
//    if (argc < 2) {
//        printf("Masa  barionowa:");
//        scanf("%d", &m_barionowa);
//    } else {
//        m_barionowa = atoi(argv[1]);
//    }

    WierszNode *lista_wierszy = NULL;

    appendWiersz(&lista_wierszy, 1.1, 67, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
    appendWiersz(&lista_wierszy, 1.25, 90, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);// wieksza
    appendWiersz(&lista_wierszy, 1.5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
    appendWiersz(&lista_wierszy, 1.12, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);// najmniejsza
    appendWiersz(&lista_wierszy, 1.3, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
    appendWiersz(&lista_wierszy, 1.11, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);

    WierszNode *przeszukane3 = findByMBar(&lista_wierszy, 1.2);

    if (przeszukane3 != NULL) printf("%f\n", przeszukane3->m_grawitacyjna);
    else {
        ParaWierszy *paraWierszy = znajdzNajblizszaPare(&lista_wierszy, 1.2);
        printf("mniejsza: %f\n", paraWierszy->mniejszy->m_grawitacyjna);
        printf("wieksza: %f\n",  paraWierszy->wiekszy->m_grawitacyjna);
    }

    destroyList(&lista_wierszy);
    return 0;
}
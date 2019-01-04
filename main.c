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
    WierszNode *wyliczony = (WierszNode *) malloc(sizeof(WierszNode));
    WierszNode *current_node = *head_ref;

    wyliczony->next = NULL;
    wyliczony->m_barionowa = m_barionowa;
    wyliczony->m_grawitacyjna = m_grawitacyjna;
    wyliczony->m_pedu = m_pedu;
    wyliczony->R = R;
    wyliczony->omega = omega;
    wyliczony->P = P;
    wyliczony->Ro_B = Ro_B;
    wyliczony->N_B = N_B;
    wyliczony->entalpia = entalpia;
    wyliczony->r_eq = r_eq;
    wyliczony->s_r = s_r;
    wyliczony->s_e = s_e;


    if (*head_ref == NULL) {
        *head_ref = wyliczony;
        return;
    }

    while (current_node->next != NULL)
        current_node = current_node->next;

    current_node->next = wyliczony;
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
    double m_bar_szukana;
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
    paraWierszy->m_bar_szukana = m_bar_szukana;

    return paraWierszy;
}

double liniowaZaleznosc(double min_y, double max_y, double min_x, double max_x) {
    double a, a1, a2, b, b1;
    a1 = (max_y - min_y);

    // wyliczanie masy grawitacyjnej
    a2 = (max_x - min_x);
    a = a1 / a2;

    b1 = ((min_y * max_x) - (max_y * min_x));
    b = b1 / a2;

    return (min_y - b) / a;
}

WierszNode *wyliczWierszLiniowo(ParaWierszy *paraWierszy) {
    WierszNode *wyliczony = malloc(sizeof(WierszNode));
    WierszNode *mniejszy = paraWierszy->mniejszy;
    WierszNode *wiekszy = paraWierszy->wiekszy;

    wyliczony->next = NULL;
    wyliczony->m_barionowa = paraWierszy->m_bar_szukana;
    wyliczony->m_grawitacyjna = liniowaZaleznosc(mniejszy->m_barionowa, wiekszy->m_barionowa, mniejszy->m_grawitacyjna, wiekszy->m_grawitacyjna);
    wyliczony->m_pedu = liniowaZaleznosc(mniejszy->m_barionowa, wiekszy->m_barionowa, mniejszy->m_pedu, wiekszy->m_pedu);
    wyliczony->R = liniowaZaleznosc(mniejszy->m_barionowa, wiekszy->m_barionowa, mniejszy->R, wiekszy->R);
    wyliczony->omega = liniowaZaleznosc(mniejszy->m_barionowa, wiekszy->m_barionowa, mniejszy->omega, wiekszy->omega);
    wyliczony->P = liniowaZaleznosc(mniejszy->m_barionowa, wiekszy->m_barionowa, mniejszy->P, wiekszy->P);
    wyliczony->Ro_B = liniowaZaleznosc(mniejszy->m_barionowa, wiekszy->m_barionowa, mniejszy->Ro_B, wiekszy->Ro_B);
    wyliczony->N_B = liniowaZaleznosc(mniejszy->m_barionowa, wiekszy->m_barionowa, mniejszy->N_B, wiekszy->N_B);
    wyliczony->entalpia = liniowaZaleznosc(mniejszy->m_barionowa, wiekszy->m_barionowa, mniejszy->entalpia, wiekszy->entalpia);
    wyliczony->r_eq = liniowaZaleznosc(mniejszy->m_barionowa, wiekszy->m_barionowa, mniejszy->r_eq, wiekszy->r_eq);
    wyliczony->s_r = liniowaZaleznosc(mniejszy->m_barionowa, wiekszy->m_barionowa, mniejszy->s_r, wiekszy->s_r);
    wyliczony->s_e = liniowaZaleznosc(mniejszy->m_barionowa, wiekszy->m_barionowa, mniejszy->s_e, wiekszy->s_e);

    return wyliczony;
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

    appendWiersz(&lista_wierszy, 1.1, 67, 23, 24, 12, 65, 34, 23, 2, 76, 23, 1);
    appendWiersz(&lista_wierszy, 1.25, 90, 24, 34, 43, 34, 23, -1, 34, 65, 54, 23);// wieksza
    appendWiersz(&lista_wierszy, 1.5, 23, 34, 56, 26, 26, 34, 14, 54, 61, 24, 34);
    appendWiersz(&lista_wierszy, 1.12, 8, 56, 65, 34, 34, 34, 45, 75, 75, 65, 43);// najmniejsza
    appendWiersz(&lista_wierszy, 1.3, 10, 86, 723, -23, 34, 23, 25, 71, 734, 82, 23);
    appendWiersz(&lista_wierszy, 1.11, 83, 25, 90, 73, 71, 26, 234, 45, 75, 75, 234);

    WierszNode *przeszukane3 = findByMBar(&lista_wierszy, 1.2);

    if (przeszukane3 != NULL) printf("%f\n", przeszukane3->m_grawitacyjna);
    else {
        ParaWierszy *paraWierszy = znajdzNajblizszaPare(&lista_wierszy, 1.2);
        printf("mniejsza: %f\n", paraWierszy->mniejszy->m_grawitacyjna);
        printf("wieksza: %f\n", paraWierszy->wiekszy->m_grawitacyjna);
    }

    destroyList(&lista_wierszy);
    return 0;
}
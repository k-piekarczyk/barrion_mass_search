//
// Created by Krzysiek on 04.01.2019.
//

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#define BUFOR_ROZNICY 99999999
#define BUFOR_LINI 256

/** Definicja listy wierszy */
typedef struct _wiersz_node {
    int generated;

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
    wyliczony->generated = 0;
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
        fprintf(stderr, "Nie znaleziono wiersza z mniejsza masa barionowa niz: %f\n", m_bar_szukana);
        exit(EXIT_FAILURE);
    } else if (wieksza == NULL) {
        fprintf(stderr, "Nie znaleziono wiersza z wieksza masa barionowa niz: %f\n", m_bar_szukana);
        exit(EXIT_FAILURE);
    }

    ParaWierszy *paraWierszy = malloc(sizeof(ParaWierszy));
    paraWierszy->mniejszy = mniejsza;
    paraWierszy->wiekszy = wieksza;
    paraWierszy->m_bar_szukana = m_bar_szukana;

    return paraWierszy;
}

double liniowaZaleznosc(double desirable_y, double min_y, double max_y, double min_x, double max_x) {

    if (min_x == max_x) return min_x;

    double a = 0, a1 = 0, a2 = 0, b = 0, b1 = 0, solution = 0;
    a1 = (max_y - min_y);

    a2 = (max_x - min_x);
    a = a1 / a2;

    b1 = ((min_y * max_x) - (max_y * min_x));
    b = b1 / a2;

    solution = (desirable_y - b) / a;
    return solution;
}

WierszNode *wyliczWierszLiniowo(ParaWierszy *paraWierszy) {
    WierszNode *wyliczony = malloc(sizeof(WierszNode));
    WierszNode *mniejszy = paraWierszy->mniejszy;
    WierszNode *wiekszy = paraWierszy->wiekszy;

    wyliczony->next = NULL;
    wyliczony->generated = 1;
    wyliczony->m_barionowa = paraWierszy->m_bar_szukana;
    wyliczony->m_grawitacyjna = liniowaZaleznosc(wyliczony->m_barionowa, mniejszy->m_barionowa, wiekszy->m_barionowa,
                                                 mniejszy->m_grawitacyjna, wiekszy->m_grawitacyjna);
    wyliczony->m_pedu = liniowaZaleznosc(wyliczony->m_barionowa, mniejszy->m_barionowa, wiekszy->m_barionowa,
                                         mniejszy->m_pedu, wiekszy->m_pedu);
    wyliczony->R = liniowaZaleznosc(wyliczony->m_barionowa, mniejszy->m_barionowa, wiekszy->m_barionowa, mniejszy->R,
                                    wiekszy->R);
    wyliczony->omega = liniowaZaleznosc(wyliczony->m_barionowa, mniejszy->m_barionowa, wiekszy->m_barionowa,
                                        mniejszy->omega, wiekszy->omega);
    wyliczony->P = liniowaZaleznosc(wyliczony->m_barionowa, mniejszy->m_barionowa, wiekszy->m_barionowa, mniejszy->P,
                                    wiekszy->P);
    wyliczony->Ro_B = liniowaZaleznosc(wyliczony->m_barionowa, mniejszy->m_barionowa, wiekszy->m_barionowa,
                                       mniejszy->Ro_B, wiekszy->Ro_B);
    wyliczony->N_B = liniowaZaleznosc(wyliczony->m_barionowa, mniejszy->m_barionowa, wiekszy->m_barionowa,
                                      mniejszy->N_B, wiekszy->N_B);
    wyliczony->entalpia = liniowaZaleznosc(wyliczony->m_barionowa, mniejszy->m_barionowa, wiekszy->m_barionowa,
                                           mniejszy->entalpia, wiekszy->entalpia);
    wyliczony->r_eq = liniowaZaleznosc(wyliczony->m_barionowa, mniejszy->m_barionowa, wiekszy->m_barionowa,
                                       mniejszy->r_eq, wiekszy->r_eq);
    wyliczony->s_r = liniowaZaleznosc(wyliczony->m_barionowa, mniejszy->m_barionowa, wiekszy->m_barionowa,
                                      mniejszy->s_r, wiekszy->s_r);
    wyliczony->s_e = liniowaZaleznosc(wyliczony->m_barionowa, mniejszy->m_barionowa, wiekszy->m_barionowa,
                                      mniejszy->s_e, wiekszy->s_e);

    return wyliczony;
}

WierszNode *findEquivalent(WierszNode **head_ref, double m_bar_szukana) {
    ParaWierszy *paraWierszy = znajdzNajblizszaPare(head_ref, m_bar_szukana);
    return wyliczWierszLiniowo(paraWierszy);
}

WierszNode *findAnswer(WierszNode **head_ref, double m_bar_szukana) {
    WierszNode *ans = findByMBar(head_ref, m_bar_szukana);
    if (ans != NULL) return ans;
    return findEquivalent(head_ref, m_bar_szukana);
}

/** Wczytywanie pliku */
void parseLine(const char *src, WierszNode **head_ref) {
    char *bufor[BUFOR_LINI];
    char *line = malloc(strlen(src) + 1);
    strcpy(line, src);

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

    char *ptr = strtok(line, ",");
    m_barionowa = atof(ptr);

    ptr = strtok(NULL, ",");
    m_grawitacyjna = atof(ptr);

    ptr = strtok(NULL, ",");
    m_pedu = atof(ptr);

    ptr = strtok(NULL, ",");
    R = atof(ptr);

    ptr = strtok(NULL, ",");
    omega = atof(ptr);

    ptr = strtok(NULL, ",");
    P = atof(ptr);

    ptr = strtok(NULL, ",");
    Ro_B = atof(ptr);

    ptr = strtok(NULL, ",");
    N_B = atof(ptr);

    ptr = strtok(NULL, ",");
    entalpia = atof(ptr);

    ptr = strtok(NULL, ",");
    r_eq = atof(ptr);

    ptr = strtok(NULL, ",");
    s_r = atof(ptr);

    ptr = strtok(NULL, ",");
    s_e = atof(ptr);

    appendWiersz(head_ref, m_barionowa, m_grawitacyjna, m_pedu, R, omega, P, Ro_B, N_B, entalpia, r_eq, s_r, s_e);
}

WierszNode *loadData(const char *fileName) {
    FILE *fp;
    char line[BUFOR_LINI];
    WierszNode *dataList = NULL;

    fp = fopen(fileName, "r");
    if (fp == NULL) {
        printf("No file '%s'.\n", fileName);
        exit(EXIT_FAILURE);
    }

    while (fgets(line, sizeof line, fp) != NULL) {
        strtok(line, "\n");
        parseLine(line, &dataList);
    }
    fclose(fp);
    return dataList;
}

/** Wypisywanie poszukiwanych danych */
void printCalculatedData(WierszNode *calcData) {
    printf("%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n", calcData->m_barionowa,
           calcData->m_grawitacyjna, calcData->m_pedu, calcData->R, calcData->omega, calcData->P, calcData->Ro_B,
           calcData->N_B, calcData->entalpia, calcData->r_eq, calcData->s_r, calcData->s_e);
}

int main(int argc, char **argv) {
    double m_barionowa;
    char *fileName;

    if (argc < 2) {
        printf("Nie podales nazwy pliku.\n");
        exit(EXIT_FAILURE);

    } else {
        fileName = argv[1];
    }

    WierszNode *lista_wierszy = loadData(fileName);

    if (argc < 3) {
        printf("Podaj mase barionowa: ");
        scanf("%f", &m_barionowa);
    } else {
        m_barionowa = atof(argv[2]);
    }

    WierszNode *answer = findAnswer(&lista_wierszy, m_barionowa);
    printCalculatedData(answer);

    if (answer->generated == 1) free(answer);
    destroyList(&lista_wierszy);
    return 0;
}
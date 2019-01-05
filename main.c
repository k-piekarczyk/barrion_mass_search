//
// Created by Krzysiek on 04.01.2019.
//

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#define DIFF_BUFFER 99999999
#define LINE_BUFFER 256

/** Definicja listy wierszy */

/**
 * Type: DataLineNode
 * Represents a line of data
 * Is a node of a one-way linked list
 * */
typedef struct _data_line_node_s {
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

    struct _data_line_node_s *next;
} DataLineNode;

/**
 * Appends a line to DataLine list
 * */
void
appendLine(DataLineNode **head_ref, double m_barionowa, double m_grawitacyjna, double m_pedu, double R, double omega,
           double P, double Ro_B, double N_B, double entalpia, double r_eq, double s_r, double s_e) {
    DataLineNode *calculated = (DataLineNode *) malloc(sizeof(DataLineNode));
    DataLineNode *current_node = *head_ref;

    calculated->next = NULL;
    calculated->generated = 0;
    calculated->m_barionowa = m_barionowa;
    calculated->m_grawitacyjna = m_grawitacyjna;
    calculated->m_pedu = m_pedu;
    calculated->R = R;
    calculated->omega = omega;
    calculated->P = P;
    calculated->Ro_B = Ro_B;
    calculated->N_B = N_B;
    calculated->entalpia = entalpia;
    calculated->r_eq = r_eq;
    calculated->s_r = s_r;
    calculated->s_e = s_e;


    if (*head_ref == NULL) {
        *head_ref = calculated;
        return;
    }

    while (current_node->next != NULL)
        current_node = current_node->next;

    current_node->next = calculated;
}

/**
 * Searches DataLine list for a Node by baryon mass
 * @returns {DataLineNode or NULL} - Node containing data or NULL
 * */
DataLineNode *findByBaryonMass(DataLineNode **head_ref, double m_bar) {
    DataLineNode *current_node = *head_ref;

    while (current_node != NULL) {
        if (current_node->m_barionowa == m_bar) return current_node;
        current_node = current_node->next;
    }

    return NULL;
}

/**
 * Frees a list memory
 * */
void destroyList(DataLineNode **head_ref) {
    DataLineNode *current_node = *head_ref;
    DataLineNode *next;

    while (current_node != NULL) {
        next = current_node->next;
        free(current_node);
        current_node = next;
    }

    *head_ref = NULL;
}

/** Przeszukiwanie najbliższych */

/**
 * Type: DataLinePair
 * Represents a pair of DataLineNodes
 * */
typedef struct _data_line_pair {
    DataLineNode *smaller;
    DataLineNode *bigger;
    double bar_m_queried;
} DataLinePair;

DataLinePair *findClosestPair(DataLineNode **head_ref, double bar_m_query) {
    DataLineNode *current_node = *head_ref;

    DataLineNode *smaller = NULL;
    DataLineNode *bigger = NULL;
    double min_diff = DIFF_BUFFER;
    double max_diff = DIFF_BUFFER;
    double diff;


    while (current_node != NULL) {
        diff = fabs(current_node->m_barionowa - bar_m_query);

        if ((current_node->m_barionowa < bar_m_query) && (diff < min_diff)) { // Szukanie mniejszej
            min_diff = diff;
            smaller = current_node;
        } else if ((current_node->m_barionowa > bar_m_query) && (diff < max_diff)) { // Szukanie wiekszej
            max_diff = diff;
            bigger = current_node;
        }

        current_node = current_node->next;
    }

    if (smaller == NULL) {
        fprintf(stderr, "Nie znaleziono wiersza z mniejsza masa barionowa niz: %f\n", bar_m_query);
        exit(EXIT_FAILURE);
    } else if (bigger == NULL) {
        fprintf(stderr, "Nie znaleziono wiersza z wieksza masa barionowa niz: %f\n", bar_m_query);
        exit(EXIT_FAILURE);
    }

    DataLinePair *linePair = malloc(sizeof(DataLinePair));
    linePair->smaller = smaller;
    linePair->bigger = bigger;
    linePair->bar_m_queried = bar_m_query;

    return linePair;
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

DataLineNode *wyliczWierszLiniowo(DataLinePair *paraWierszy) {
    DataLineNode *wyliczony = malloc(sizeof(DataLineNode));
    DataLineNode *mniejszy = paraWierszy->smaller;
    DataLineNode *wiekszy = paraWierszy->bigger;

    wyliczony->next = NULL;
    wyliczony->generated = 1;
    wyliczony->m_barionowa = paraWierszy->bar_m_queried;
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

DataLineNode *findEquivalent(DataLineNode **head_ref, double m_bar_szukana) {
    DataLinePair *paraWierszy = findClosestPair(head_ref, m_bar_szukana);
    return wyliczWierszLiniowo(paraWierszy);
}

DataLineNode *findAnswer(DataLineNode **head_ref, double m_bar_szukana) {
    DataLineNode *ans = findByBaryonMass(head_ref, m_bar_szukana);
    if (ans != NULL) return ans;
    return findEquivalent(head_ref, m_bar_szukana);
}

/** Wczytywanie pliku */
void parseLine(const char *src, DataLineNode **head_ref) {
    char *bufor[LINE_BUFFER];
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

    appendLine(head_ref, m_barionowa, m_grawitacyjna, m_pedu, R, omega, P, Ro_B, N_B, entalpia, r_eq, s_r, s_e);
}

DataLineNode *loadData(const char *fileName) {
    FILE *fp;
    char line[LINE_BUFFER];
    DataLineNode *dataList = NULL;

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
void printCalculatedData(DataLineNode *calcData) {
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

    DataLineNode *lista_wierszy = loadData(fileName);

    if (argc < 3) {
        printf("Podaj mase barionowa: ");
        scanf("%f", &m_barionowa);
    } else {
        m_barionowa = atof(argv[2]);
    }

    DataLineNode *answer = findAnswer(&lista_wierszy, m_barionowa);
    printCalculatedData(answer);

    if (answer->generated == 1) free(answer);
    destroyList(&lista_wierszy);
    return 0;
}
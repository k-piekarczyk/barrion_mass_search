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
    int calculated;         // informuje czy był w plikach (0) czy został wyliczony (1)

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
    calculated->calculated = 0;
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

/**
 * Finds the closest bigger-smaller pair for queried baryon mass
 * @sidefects Finishes program with EXIT_FAILIURE in case no smaller-bigger pair is not found in provided data
 * @returns {DataLinePair}
 * */
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

/**
 * Calculates a value for a certain 'y' based on a linear relation
 * @returns {double} - calculated value
 * */
double linearRelation(double desired_y, double min_y, double max_y, double min_x, double max_x) {

    if (min_x == max_x) return min_x;

    double a = 0, a1 = 0, a2 = 0, b = 0, b1 = 0, solution = 0;
    a1 = (max_y - min_y);

    a2 = (max_x - min_x);
    a = a1 / a2;

    b1 = ((min_y * max_x) - (max_y * min_x));
    b = b1 / a2;

    solution = (desired_y - b) / a;
    return solution;
}

/**
 * Calculates all values of a DataLine for a certain baryon mass based on a linear relation of a bigger-smaller pair provided
 * @returns {DataLineNode} - calculated DataLine
 * */
DataLineNode *calculateDataLineLinearly(DataLinePair *linePair) {
    DataLineNode *calculated = malloc(sizeof(DataLineNode));
    DataLineNode *smaller = linePair->smaller;
    DataLineNode *bigger = linePair->bigger;

    calculated->next = NULL;
    calculated->calculated = 1;
    calculated->m_barionowa = linePair->bar_m_queried;
    calculated->m_grawitacyjna = linearRelation(calculated->m_barionowa, smaller->m_barionowa, bigger->m_barionowa,
                                               smaller->m_grawitacyjna, bigger->m_grawitacyjna);
    calculated->m_pedu = linearRelation(calculated->m_barionowa, smaller->m_barionowa, bigger->m_barionowa,
                                       smaller->m_pedu, bigger->m_pedu);
    calculated->R = linearRelation(calculated->m_barionowa, smaller->m_barionowa, bigger->m_barionowa, smaller->R,
                                  bigger->R);
    calculated->omega = linearRelation(calculated->m_barionowa, smaller->m_barionowa, bigger->m_barionowa,
                                      smaller->omega, bigger->omega);
    calculated->P = linearRelation(calculated->m_barionowa, smaller->m_barionowa, bigger->m_barionowa, smaller->P,
                                  bigger->P);
    calculated->Ro_B = linearRelation(calculated->m_barionowa, smaller->m_barionowa, bigger->m_barionowa,
                                     smaller->Ro_B, bigger->Ro_B);
    calculated->N_B = linearRelation(calculated->m_barionowa, smaller->m_barionowa, bigger->m_barionowa,
                                    smaller->N_B, bigger->N_B);
    calculated->entalpia = linearRelation(calculated->m_barionowa, smaller->m_barionowa, bigger->m_barionowa,
                                         smaller->entalpia, bigger->entalpia);
    calculated->r_eq = linearRelation(calculated->m_barionowa, smaller->m_barionowa, bigger->m_barionowa,
                                     smaller->r_eq, bigger->r_eq);
    calculated->s_r = linearRelation(calculated->m_barionowa, smaller->m_barionowa, bigger->m_barionowa,
                                    smaller->s_r, bigger->s_r);
    calculated->s_e = linearRelation(calculated->m_barionowa, smaller->m_barionowa, bigger->m_barionowa,
                                    smaller->s_e, bigger->s_e);

    return calculated;
}

/**
 * Calculates a DataLine based on provided baryon mass not found in data
 * @returns {DataLineNode} - calculated DataLine
 * */
DataLineNode *findEquivalent(DataLineNode **head_ref, double bar_m_query) {
    DataLinePair *paraWierszy = findClosestPair(head_ref, bar_m_query);
    return calculateDataLineLinearly(paraWierszy);
}

/**
 * Finds or calculates a DataLine by baryon mass and provided data
 * @returns {DataLineNode}
 * */
DataLineNode *findAnswer(DataLineNode **head_ref, double bar_m_query) {
    DataLineNode *ans = findByBaryonMass(head_ref, bar_m_query);
    if (ans != NULL) return ans;
    return findEquivalent(head_ref, bar_m_query);
}

/** Wczytywanie pliku */

/**
 * Parses a string of comma separated floating-point values and appends the data to provided DataLine list
 * */
void parseLine(const char *src, DataLineNode **head_ref) {
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

/**
 * Loads provided data from file to an in memory list
 * @returns {DataLineNode} - loaded DataLine list
 * */
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

/**
 * Prints provided DataLineNode as a line of comma separated floating-point values
 * */
void printCalculatedData(DataLineNode *calcData) {
    printf("%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n", calcData->m_barionowa,
           calcData->m_grawitacyjna, calcData->m_pedu, calcData->R, calcData->omega, calcData->P, calcData->Ro_B,
           calcData->N_B, calcData->entalpia, calcData->r_eq, calcData->s_r, calcData->s_e);
}

int main(int argc, char **argv) {
    double baryon_mass;
    char *fileName;

    if (argc < 2) {
        printf("Nie podales nazwy pliku.\n");
        exit(EXIT_FAILURE);

    } else {
        fileName = argv[1];
    }

    DataLineNode *line_list = loadData(fileName);

    if (argc < 3) {
        printf("Podaj mase barionowa: ");
        scanf("%f", &baryon_mass);
    } else {
        baryon_mass = atof(argv[2]);
    }

    DataLineNode *answer = findAnswer(&line_list, baryon_mass);
    printCalculatedData(answer);

    if (answer->calculated == 1) free(answer);
    destroyList(&line_list);
    return 0;
}
/****************************************************************************************************
* Region growing tutorial: http://notmatthancock.github.io/2017/10/09/region-growing-wrapping-c.html
****************************************************************************************************/
#ifndef STACK_H
#define STACK_H

#include <stdlib.h>
#include <stdio.h>

//Element type
typedef struct element {
    struct element * next;
    int i,j,k;
} element;

//Element stack
typedef struct stack {
    int n_elements;
    element * top;
} stack;


/*****************************************
* Initialize the stack
*****************************************/
void stack_init(stack * S) {
    S->n_elements = 0;
    S->top = NULL;
}

/*****************************************
* Initialize the element
*****************************************/
void element_init(element * el, int i, int j, int k) {
    el->next = NULL;
    el->i = i;
    el->j = j;
    el->k = k;
}

/*****************************************
* Push an element to the stack structure
*****************************************/
void stack_push(stack * S, element * el) {
    S->n_elements++;
    el->next = (S->top);
    S->top = el;
}

/*****************************************
* Pop an element from the stack structure
*****************************************/
element * stack_pop(stack * S) {
    element * el;
    if (S->n_elements == 0) {
        return NULL;
    }
    else {
        S->n_elements--;
        el = S->top;
        if (S->n_elements == 0)
            S->top = NULL;
        else
            S->top = el->next;

        return el;
    }
}
#endif
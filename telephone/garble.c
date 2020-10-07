#include "garble.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


/*
 * GARBLE simulates the mispronunciation and blending of syllables
 * than young children might possess.  Feel free to modify garble as
 * you see fit!
 */
void garble(char* message)
{
    int sflag = rand() % 4;
    switch (sflag) {

        /* Add a random lower case letter */
        case 0:
        {
            int n = strlen(message);
            int k = rand() % n;
            message[k] = 'a' + (rand() % 26);
            break;
        }

        /* Add a random upper case letter */
        case 1:
        {
            int n = strlen(message);
            int k = rand() % n;
            message[k] = 'A' + (rand() % 26);
            break;
        }

        /* Add a random "blegh" */
        case 2:
        {
            int n = strlen(message);
            if (n > 5)
                memcpy(message + (rand() % (n-5)), "blegh", 5);
            break;
        }

        /* Add a random "gargh" */
        case 3:
        {
            int n = strlen(message);
            if (n > 5)
                memcpy(message + (rand() % (n-5)), "gargh", 5);
            break;
        }
    }
}

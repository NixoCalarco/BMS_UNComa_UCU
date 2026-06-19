#include "freertos/idf_additions.h"
#include "freertos/projdefs.h"
#include <esp_sleep.h>
#include <freertos/FreeRTOS.h>
#include <stdio.h>

#include "kalman.h"

// static QueueHandle_t msg_queue;
// static const uint8_t msg_queue_len = 40;

// static void consumer(void *arg) {
//     KalmanVals values;
//
//     int to_wait_ms = 1000; // the maximal blocking waiting time of millisecond
//     const TickType_t xTicksToWait = pdMS_TO_TICKS(to_wait_ms);
//
//     while (true) {
//         if (xQueueReceive(msg_queue, &values, xTicksToWait))
//             printf("%s, "
//                    "V, %f "
//                    "Se, %f "
//                    "XXe, %f "
//                    "TT, %f "
//                    "TE, %f\n",
//                    values.algorithm == EXTENDED ? "Kalman Extended" : "Kalman",
//                    values.V,
//                    values.Se,
//                    values.XXe,
//                    values.TT,
//                    values.TE);
//         else
//             printf("Esperando datos...\n");
//
//         vTaskDelay(pdMS_TO_TICKS(1000));
//     }
//
//     vTaskDelete(nullptr);
// }

void push_to_queue(KalmanVals values) {
    // const TickType_t xTicksToWait = pdMS_TO_TICKS(1000);
    // if (!xQueueSend(msg_queue, &values, xTicksToWait))
    //     printf("Queue full!!\n");

    // delay acá para que tranque kalman()
    if (esp_sleep_enable_timer_wakeup(1'000'000) != ESP_OK)
        printf("No se pudo setear timer.\n");
    if (esp_light_sleep_start() != ESP_OK)
        printf("No puedo dormir\n");

    // vTaskDelay(pdMS_TO_TICKS(1000));
}

static void producer(void *arg) {
    while (true) {
        kalman(push_to_queue);
    // printf("Hola!\n");
    }

    vTaskDelete(nullptr);
}

void app_main(void) {
    printf("Hola!\n");

    // msg_queue = xQueueCreate(msg_queue_len, sizeof(KalmanVals));

    // if (esp_sleep_enable_timer_wakeup(1'000'000) != ESP_OK)
    //     printf("No se pudo setear timer.\n");

    // xTaskCreate(consumer, "consumer", 4096, nullptr, 3, nullptr);
    xTaskCreate(producer, "producer", 4096, nullptr, 3, nullptr);
}

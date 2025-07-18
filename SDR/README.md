# DVB-S2 Transceiver Using SDR
## Experiment Procedure
1. **Set the Radio IDs**
    
    For simplicity, make sure that the radio devices are correctly assigned as follow:
    
    ```matlab
    radioTx.RadioID = "usb:1"; % Right side (Transmitter)
    radioRx.RadioID = "usb:0"; % Left side (Receiver)
    ```
    
2. **Start the Receiver first**
To ensure the correct radio device is captured, run the **receiver program first**, followed by the **transmitter**.
3. **Transmit and Capture**
    - The transmitter continuously sends the same signal 1000 times.
    - Start the receiver at any time during this process to capture the transmission.
## Reference
MATLAB Documentations:
1. [DVB-S2 Transmitter Using Software-Defined Radio](https://www.mathworks.com/help/satcom/ug/dvbs2_transmitter_using_sdr.html)
2. [DVB-S2 Receiver Using Software-Defined Radio](https://www.mathworks.com/help/satcom/ug/dvbs2_receiver_using_sdr.html)

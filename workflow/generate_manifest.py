## Edit these values with your min and max sample number and sample intials
## Current values for this project
MIN_SAMPLE_NUM = 1
MAX_SAMPLE_NUM = 61
POOL_INITIALS = "KD"

## Do not edit SWITCH or RESET_VALUE
SWITCH = 0
RESET_VALUE = 0

RESET_VALUE = MIN_SAMPLE_NUM

print('SampleID')
while SWITCH == 0:
    if MIN_SAMPLE_NUM < 10:
        print(POOL_INITIALS +"00"+ str(MIN_SAMPLE_NUM))
        print(POOL_INITIALS +"00"+ str(MIN_SAMPLE_NUM))

    elif MIN_SAMPLE_NUM >= 10:
        print(POOL_INITIALS + "0" + str(MIN_SAMPLE_NUM))
        print(POOL_INITIALS + "0" + str(MIN_SAMPLE_NUM))

    MIN_SAMPLE_NUM += 1
    if MIN_SAMPLE_NUM > MAX_SAMPLE_NUM:
        SWITCH = 1

MIN_SAMPLE_NUM = RESET_VALUE
SWITCH = 0
print('\n\nabsolute-file-path')
while SWITCH == 0:
    if MIN_SAMPLE_NUM < 10:
        print("$PWD/demux_seqs/" + POOL_INITIALS +"00"+ str(MIN_SAMPLE_NUM) + "_R1.fastq.gz")
        print("$PWD/demux_seqs/" + POOL_INITIALS +"00"+ str(MIN_SAMPLE_NUM) + "_R2.fastq.gz")

    elif MIN_SAMPLE_NUM >= 10:
        print("$PWD/demux_seqs/" + POOL_INITIALS +"0"+ str(MIN_SAMPLE_NUM) + "_R1.fastq.gz")
        print("$PWD/demux_seqs/" + POOL_INITIALS +"0"+ str(MIN_SAMPLE_NUM) + "_R2.fastq.gz")

    MIN_SAMPLE_NUM += 1
    if MIN_SAMPLE_NUM > MAX_SAMPLE_NUM:
        SWITCH = 1

MIN_SAMPLE_NUM = RESET_VALUE
SWITCH = 0
print('\ndirection')
while SWITCH == 0:
    if MIN_SAMPLE_NUM <= MAX_SAMPLE_NUM:
        print("forward")
        print("reverse")
        MIN_SAMPLE_NUM += 1

    else:
        SWITCH = 1


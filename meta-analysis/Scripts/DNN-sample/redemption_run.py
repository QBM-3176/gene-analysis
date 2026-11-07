import os
from dotenv import load_dotenv
from deep_neural_network import DeepNeuralNetworkRunner

# environment

LITE_MODE = False

load_dotenv(override=True)
hf_token = os.environ['HF_TOKEN']
login(hf_token, add_to_git_credential=True)

username = "ed-donner"
dataset = f"{username}/items_lite" if LITE_MODE else f"{username}/items_full"

train, val, test = Item.from_hub(dataset)

print(f"Loaded {len(train):,} training items, {len(val):,} validation items, {len(test):,} test items")

runner = DeepNeuralNetworkRunner(train, val[:1000])
runner.setup()

#runner.load('deep_neural_network.pth')                          # to load on a MAC
runner.load('deep_neural_network.pth', 'cpu')                   # use this one if want to load on a Win PC device --- if so comment out the line above

def deep_neural_network(item):
    return runner.inference(item)

evaluate(deep_neural_network, test)

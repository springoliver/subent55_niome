module.exports = {
  apps: [
    // {
    //   name: "springhot",
    //   script: "python",
    //   args: [
    //     "-m", "neurons.miner",
    //     "--netuid", "55",
    //     "--subtensor.network", "finney",
    //     "--wallet.name", "spring",
    //     "--wallet.hotkey", "springhot",
    //     "--axon.port", "8105",
    //     "--blacklist.force_validator_permit",
    //     "--logging.debug"
    //   ],
    //   cwd: ".",
    //   interpreter: "none",
    //   autorestart: true,
    //   restart_delay: 5000
    // },
    // {
    //   name: "spring01",
    //   script: "python",
    //   args: [
    //     "-m", "neurons.miner",
    //     "--netuid", "55",
    //     "--subtensor.network", "finney",
    //     "--wallet.name", "spring",
    //     "--wallet.hotkey", "spring01",
    //     "--axon.port", "8101",
    //     "--logging.debug"
    //   ],
    //   cwd: ".",
    //   interpreter: "none",
    //   autorestart: true,
    //   restart_delay: 5000
    // },
    // {
    //   name: "spring02",
    //   script: "python",
    //   args: [
    //     "-m", "neurons.miner",
    //     "--netuid", "55",
    //     "--subtensor.network", "finney",
    //     "--wallet.name", "spring",
    //     "--wallet.hotkey", "spring02",
    //     "--axon.port", "8112",
    //     "--blacklist.force_validator_permit",
    //     "--logging.debug"
    //   ],
    //   cwd: ".",
    //   interpreter: "none",
    //   autorestart: true,
    //   restart_delay: 5000
    // },
    {
      name: "spring03",
      script: "python",
      args: [
        "-m", "neurons.miner",
        "--netuid", "55",
        "--subtensor.network", "finney",
        "--wallet.name", "spring",
        "--wallet.hotkey", "spring03",
        "--axon.port", "8103",
        "--logging.debug"
      ],
      cwd: ".",
      interpreter: "none",
      autorestart: true,
      restart_delay: 5000
    },
    {
      name: "spring04",
      script: "python",
      args: [
        "-m", "neurons.miner",
        "--netuid", "55",
        "--subtensor.network", "finney",
        "--wallet.name", "spring",
        "--wallet.hotkey", "spring04",
        "--axon.port", "8104",
        "--blacklist.force_validator_permit",
        "--logging.debug"
      ],
      cwd: ".",
      interpreter: "none",
      autorestart: true,
      restart_delay: 5000
    },
    {
      name: "spring05",
      script: "python",
      args: [
        "-m", "neurons.miner",
        "--netuid", "55",
        "--subtensor.network", "finney",
        "--wallet.name", "spring",
        "--wallet.hotkey", "spring05",
        "--axon.port", "8107",
        "--logging.debug"
      ],
      cwd: ".",
      interpreter: "none",
      autorestart: true,
      restart_delay: 5000
    },
    // {
    //   name: "spring06",
    //   script: "python",
    //   args: [
    //     "-m", "neurons.miner",
    //     "--netuid", "55",
    //     "--subtensor.network", "finney",
    //     "--wallet.name", "spring",
    //     "--wallet.hotkey", "spring06",
    //     "--axon.port", "8110",
    //     "--blacklist.force_validator_permit",
    //     "--logging.debug"
    //   ],
    //   cwd: ".",
    //   interpreter: "none",
    //   autorestart: true,
    //   restart_delay: 5000
    // },
    // {
    //   name: "spring07",
    //   script: "python",
    //   args: [
    //     "-m", "neurons.miner",
    //     "--netuid", "55",
    //     "--subtensor.network", "finney",
    //     "--wallet.name", "spring",
    //     "--wallet.hotkey", "spring07",
    //     "--axon.port", "8111",
    //     "--logging.debug"
    //   ],
    //   cwd: ".",
    //   interpreter: "none",
    //   autorestart: true,
    //   restart_delay: 5000
    // },
    // {
    //   name: "spring08",
    //   script: "python",
    //   args: [
    //     "-m", "neurons.miner",
    //     "--netuid", "55",
    //     "--subtensor.network", "finney",
    //     "--wallet.name", "spring",
    //     "--wallet.hotkey", "spring08",
    //     "--axon.port", "8108",
    //     "--blacklist.force_validator_permit",
    //     "--logging.debug"
    //   ],
    //   cwd: ".",
    //   interpreter: "none",
    //   autorestart: true,
    //   restart_delay: 5000
    // },
    // {
    //   name: "spring09",
    //   script: "python",
    //   args: [
    //     "-m", "neurons.miner",
    //     "--netuid", "55",
    //     "--subtensor.network", "finney",
    //     "--wallet.name", "spring",
    //     "--wallet.hotkey", "spring09",
    //     "--axon.port", "8109",
    //     "--logging.debug"
    //   ],
    //   cwd: ".",
    //   interpreter: "none",
    //   autorestart: true,
    //   restart_delay: 5000
    // },
    {
      name: "dragon",
      script: "python",
      args: [
        "-m", "neurons.miner",
        "--netuid", "55",
        "--subtensor.network", "finney",
        "--wallet.name", "dragon",
        "--wallet.hotkey", "dragon",
        "--axon.port", "8106",
        "--blacklist.force_validator_permit",
        "--logging.debug"
      ],
      cwd: ".",
      interpreter: "none",
      autorestart: true,
      restart_delay: 5000
    }
  ]
};

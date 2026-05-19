module.exports = {
  apps: [
    {
      name: "miner-1",
      script: "python",
      args: [
        "-m", "neurons.miner",
        "--netuid", "55",
        "--subtensor.network", "finney",
        "--wallet.name", "pang8512",
        "--wallet.hotkey", "miner-pang",
        "--axon.port", "8091",
        "--logging.debug"
      ],
      cwd: ".",
      interpreter: "none",
      autorestart: true,
      restart_delay: 5000
    }
  ]
};

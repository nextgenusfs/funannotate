#!/usr/bin/env bash
# Utility to toggle Rust PASA tools on/off for ad-hoc testing
#
# Mirrors toggle_trinity_rust.sh's mechanism (rename-to-.disabled, reversible)
# for PASA's Rust binaries, which live under $CONDA_PREFIX/opt/pasa/bin (see
# install_scripts/pixi_install_pasa.sh), not $CONDA_PREFIX/bin like Trinity's.
#
# Usage:
#   toggle_pasa_rust.sh status              # Show current status
#   toggle_pasa_rust.sh disable             # Hide Rust tools (use standard PASA pipeline)
#   toggle_pasa_rust.sh enable              # Show Rust tools (use Rust)
#   toggle_pasa_rust.sh reset               # Restore to original state

set -euo pipefail

# Check if in a pixi/conda environment
if [ -z "${CONDA_PREFIX:-}" ]; then
    echo "ERROR: Not in a pixi/conda environment"
    echo "Please run: pixi shell"
    exit 1
fi

BIN_DIR="${CONDA_PREFIX}/opt/pasa/bin"
DISABLE_SUFFIX=".disabled"

# Rust PASA tools (see install_scripts/pixi_install_pasa.sh's own completeness
# check for this same list of names)
PASA_TOOLS=(
    "pasa_rust"
    "slclust_rust"
    "cdbyank_rust"
    "faidx_rust"
)

print_status() {
    echo "Rust PASA Tools Status"
    echo "========================"
    echo "CONDA_PREFIX: $CONDA_PREFIX"
    echo "BIN_DIR: $BIN_DIR"
    echo ""
    echo "Tools:"
    for tool in "${PASA_TOOLS[@]}"; do
        if [ -x "$BIN_DIR/$tool" ]; then
            echo "  ✓ $tool (enabled)"
        elif [ -f "$BIN_DIR/${tool}${DISABLE_SUFFIX}" ]; then
            echo "  ✗ $tool (disabled - using standard PASA pipeline)"
        else
            echo "  ? $tool (not found)"
        fi
    done
    echo ""
}

disable_rust() {
    echo "Disabling Rust PASA tools (will use standard PASA pipeline)..."
    local disabled_count=0

    for tool in "${PASA_TOOLS[@]}"; do
        if [ -x "$BIN_DIR/$tool" ]; then
            mv "$BIN_DIR/$tool" "$BIN_DIR/${tool}${DISABLE_SUFFIX}"
            echo "  Disabled: $tool"
            disabled_count=$((disabled_count + 1))
        elif [ -f "$BIN_DIR/${tool}${DISABLE_SUFFIX}" ]; then
            echo "  Already disabled: $tool"
        fi
    done

    if [ $disabled_count -gt 0 ]; then
        echo ""
        echo "✓ Disabled $disabled_count Rust tool(s)"
        echo "PASA will now use its standard pipeline"
        echo ""
        echo "To re-enable later, run: toggle_pasa_rust.sh enable"
    fi
}

enable_rust() {
    echo "Enabling Rust PASA tools..."
    local enabled_count=0

    for tool in "${PASA_TOOLS[@]}"; do
        if [ -f "$BIN_DIR/${tool}${DISABLE_SUFFIX}" ]; then
            mv "$BIN_DIR/${tool}${DISABLE_SUFFIX}" "$BIN_DIR/$tool"
            chmod +x "$BIN_DIR/$tool"
            echo "  Enabled: $tool"
            enabled_count=$((enabled_count + 1))
        elif [ -x "$BIN_DIR/$tool" ]; then
            echo "  Already enabled: $tool"
        fi
    done

    if [ $enabled_count -gt 0 ]; then
        echo ""
        echo "✓ Enabled $enabled_count Rust tool(s)"
        echo "PASA will now use optimized Rust versions"
        echo ""
        echo "To disable later, run: toggle_pasa_rust.sh disable"
    fi
}

reset_state() {
    echo "Resetting to original state..."
    enable_rust  # Re-enable any disabled tools
}

# Main
case "${1:-status}" in
    status)
        print_status
        ;;
    disable)
        disable_rust
        print_status
        ;;
    enable)
        enable_rust
        print_status
        ;;
    reset)
        reset_state
        print_status
        ;;
    *)
        echo "Usage: $0 {status|enable|disable|reset}"
        echo ""
        echo "Commands:"
        echo "  status   - Show current status of Rust PASA tools"
        echo "  enable   - Enable Rust PASA tools (use optimized versions)"
        echo "  disable  - Disable Rust PASA tools (use standard pipeline)"
        echo "  reset    - Reset to default state (re-enable all tools)"
        echo ""
        echo "Examples:"
        echo "  $0 status          # Check which tools are active"
        echo "  $0 disable         # Use standard pipeline for benchmarking"
        echo "  $0 enable          # Switch back to Rust versions"
        exit 1
        ;;
esac

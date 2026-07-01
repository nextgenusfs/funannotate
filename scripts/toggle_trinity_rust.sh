#!/usr/bin/env bash
# Utility to toggle Rust Trinity tools on/off for ad-hoc testing
#
# Usage:
#   toggle_trinity_rust.sh status              # Show current status
#   toggle_trinity_rust.sh disable             # Hide Rust tools (use Perl)
#   toggle_trinity_rust.sh enable              # Show Rust tools (use Rust)
#   toggle_trinity_rust.sh reset               # Restore to original state

set -euo pipefail

# Check if in a pixi/conda environment
if [ -z "${CONDA_PREFIX:-}" ]; then
    echo "ERROR: Not in a pixi/conda environment"
    echo "Please run: pixi shell"
    exit 1
fi

BIN_DIR="${CONDA_PREFIX}/bin"
DISABLE_SUFFIX=".disabled"

# Rust Trinity tools
TRINITY_TOOLS=(
    "sam_to_read_coords"
    "extract_reads_per_partition"
    "fragment_coverage_writer"
    "define_coverage_partitions"
)

print_status() {
    echo "Rust Trinity Tools Status"
    echo "=========================="
    echo "CONDA_PREFIX: $CONDA_PREFIX"
    echo ""
    echo "Tools:"
    for tool in "${TRINITY_TOOLS[@]}"; do
        if [ -x "$BIN_DIR/$tool" ]; then
            echo "  ✓ $tool (enabled)"
        elif [ -f "$BIN_DIR/${tool}${DISABLE_SUFFIX}" ]; then
            echo "  ✗ $tool (disabled - using Perl version)"
        else
            echo "  ? $tool (not found)"
        fi
    done
    echo ""
}

disable_rust() {
    echo "Disabling Rust Trinity tools (will use Perl versions)..."
    local disabled_count=0

    for tool in "${TRINITY_TOOLS[@]}"; do
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
        echo "Trinity will now use Perl versions from conda"
        echo ""
        echo "To re-enable later, run: toggle_trinity_rust.sh enable"
    fi
}

enable_rust() {
    echo "Enabling Rust Trinity tools..."
    local enabled_count=0

    for tool in "${TRINITY_TOOLS[@]}"; do
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
        echo "Trinity will now use optimized Rust versions"
        echo ""
        echo "To disable later, run: toggle_trinity_rust.sh disable"
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
        echo "  status   - Show current status of Rust Trinity tools"
        echo "  enable   - Enable Rust Trinity tools (use optimized versions)"
        echo "  disable  - Disable Rust Trinity tools (use Perl versions)"
        echo "  reset    - Reset to default state (re-enable all tools)"
        echo ""
        echo "Examples:"
        echo "  $0 status          # Check which tools are active"
        echo "  $0 disable         # Use Perl versions for benchmarking"
        echo "  $0 enable          # Switch back to Rust versions"
        exit 1
        ;;
esac
